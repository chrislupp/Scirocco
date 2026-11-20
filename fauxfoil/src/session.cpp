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
#include "session.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "bl_profiles.h"
#include "format.h"
#include "postscript.h"
#include "scirocco/core/boundary_layer.h"
#include "scirocco/core/geometry.h"
#include "scirocco/core/paneling.h"
#include "scirocco/data/state.h"
#include "scirocco/properties/bl_properties.h"
#include "scirocco/scirocco.h"
#include "solution.h"

namespace fauxfoil {

using scirocco::FoilState;
using scirocco::Vector2d;
using scirocco::Vector4d;

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr const char* kProgram = "fauXFOIL";

/// Thrown by the prompts when input ends; unwinds to Run().
struct EndOfInput {};

/// A command line split as XFOIL does: the keyword (upper case, at most four
/// characters) and the rest of the line.
struct Command {
    std::string name;
    std::string arg;
    std::vector<double> numbers;
};

std::string Trim(const std::string& text) {
    const std::size_t begin = text.find_first_not_of(" \t");
    if (begin == std::string::npos) return "";
    const std::size_t end = text.find_last_not_of(" \t");
    return text.substr(begin, end - begin + 1);
}

std::string Upper(std::string text) {
    for (char& c : text)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return text;
}

/// Leading numbers of an argument string (blank or comma separated).
std::vector<double> Numbers(const std::string& arg) {
    std::vector<double> values;
    std::string text = arg;
    std::replace(text.begin(), text.end(), ',', ' ');
    std::istringstream stream(text);
    std::string token;
    while (stream >> token) {
        char* end = nullptr;
        const double value = std::strtod(token.c_str(), &end);
        if (end == token.c_str() || *end != '\0' || !std::isfinite(value))
            break;
        values.push_back(value);
    }
    return values;
}

Command Parse(const std::string& line) {
    Command command;
    const std::string text = Trim(line);
    const std::size_t split = text.find_first_of(" \t");
    command.name = Upper(text.substr(0, split)).substr(0, 4);
    if (split != std::string::npos) command.arg = Trim(text.substr(split));
    command.numbers = Numbers(command.arg);
    return command;
}

bool Is(const Command& command, std::initializer_list<const char*> names) {
    for (const char* name : names)
        if (command.name == name) return true;
    return false;
}

bool FileExists(const std::string& path) {
    std::ifstream stream(path);
    return static_cast<bool>(stream);
}

/// Freestream turbulence level [%] for an Ncrit, by Mack's correlation.
double TurbulenceLevel(double ncrit) {
    return 100.0 * std::exp(-(ncrit + 8.43) / 2.4);
}

/// x/c of transition on one surface and where it lies among the nodes.
struct Transition {
    bool free = false; /**< Free transition (else the trailing edge) */
    double x = 1.0;
    int station = 0;   /**< XFOIL boundary-layer station index */
    double node = 0.0; /**< Fractional node index counted from the
                            upper trailing edge, as XFOIL's Itr */
};

/// Transition on surface @p si (0 lower, 1 upper). The location within the
/// transition interval solves the same amplification equation as
/// ResidualTransition: ncrit - n1 - (xt - x1) * (damp1 + damp_t) / 2 = 0.
Transition FindTransition(const FoilState& state, int si) {
    Transition result;
    const auto& stations = state.vsol.Is[si];
    const int n_nodes = state.foil.N;
    if (stations.empty()) return result;
    std::size_t k = 0;
    while (k < stations.size() &&
           !(static_cast<std::size_t>(stations[k]) < state.vsol.turb.size() &&
             state.vsol.turb[stations[k]]))
        ++k;
    if (k == stations.size()) {
        const int te = stations.back();
        result.x = state.foil.x(0, te);
        result.station = static_cast<int>(stations.size()) + 1;
        result.node = n_nodes - te;
        return result;
    }
    result.free = true;
    result.station = static_cast<int>(k) + 2;
    if (k == 0) {
        result.x = state.foil.x(0, stations[0]);
        result.node = n_nodes - stations[0];
        return result;
    }
    const int i1 = stations[k - 1];
    const int i2 = stations[k];
    const double x1 = state.isol.xi(i1);
    const double x2 = state.isol.xi(i2);
    const double dx = x2 - x1;
    const Vector4d U1 = state.glob.U.col(i1);
    const Vector4d U2 = state.glob.U.col(i2);
    scirocco::Param param = scirocco::BuildParam(state, si);
    scirocco::StationParam(state, param, i1);
    param.turb = false;
    const double ncrit = param.ncrit;
    double xt = x1 + 0.5 * dx;
    if (std::fabs(dx) > 1e-14) {
        constexpr int kNewton = 30;
        for (int it = 0; it < kNewton; ++it) {
            const double w2 = (xt - x1) / dx;
            Vector4d Ut = (1.0 - w2) * U1 + w2 * U2;
            Vector4d Ut_xt = (U2 - U1) / dx;
            Ut(2) = ncrit;
            Ut_xt(2) = 0.0;
            auto [damp1, damp1_U1] = scirocco::GetDamp(U1, param);
            auto [dampt, dampt_Ut] = scirocco::GetDamp(Ut, param);
            dampt_Ut(2) = 0.0;
            const double R = ncrit - U1(2) - 0.5 * (xt - x1) * (damp1 + dampt);
            if (std::fabs(R) < 1e-14) break;
            const double R_xt =
                -0.5 * (damp1 + dampt) - 0.5 * (xt - x1) * dampt_Ut.dot(Ut_xt);
            if (R_xt == 0.0) break;
            double step = -R / R_xt;
            const double dmax =
                0.2 * dx * (1.1 - static_cast<double>(it) / kNewton);
            if (std::fabs(step) > dmax) step *= dmax / std::fabs(step);
            xt += step;
        }
    }
    const double w =
        std::fabs(dx) > 1e-14 ? std::clamp((xt - x1) / dx, 0.0, 1.0) : 0.0;
    result.x = (1.0 - w) * state.foil.x(0, i1) + w * state.foil.x(0, i2);
    result.node = n_nodes - ((1.0 - w) * i1 + w * i2);
    return result;
}

const char* const kTopMenuHead = R"(
   QUIT    Exit program

  .OPER    Direct operating point(s)
  .MDES    Complex mapping design routine
  .QDES    Surface speed design routine
  .GDES    Geometry design routine

   SAVE f  Write airfoil to labeled coordinate file
   PSAV f  Write airfoil to plain coordinate file
   ISAV f  Write airfoil to ISES coordinate file
   MSAV f  Write airfoil to MSES coordinate file
   REVE    Reverse written-airfoil node ordering
   DELI i  Change written-airfoil file delimiters

   LOAD f  Read buffer airfoil from coordinate file
   NACA i  Set NACA 4,5-digit airfoil and buffer airfoil
   INTE    Set buffer airfoil by interpolating two airfoils
   NORM    Buffer airfoil normalization toggle
   HALF    Halve the number of points in buffer airfoil
)";

const char* const kTopMenuTail = R"(
  .PPAR    Show/change paneling

  .PLOP    Plotting options

   WDEF f  Write  current-settings file
   RDEF f  Reread current-settings file
   NAME s  Specify new airfoil name
   NINC    Increment name version number

   Z       Zoom    | (available in all menus)
)"
                                 "   U       Unzoom  | \n";

const char* const kOperMenu = R"(
   <cr>     Return to Top Level
   !        Redo last ALFA,CLI,CL,ASEQ,CSEQ,VELS

   Visc r   Toggle Inviscid/Viscous mode
  .VPAR     Change BL parameter(s)
   Re   r   Change Reynolds number
   Mach r   Change Mach number
   Type i   Change type of Mach,Re variation with CL
   ITER     Change viscous-solution iteration limit
   INIT     Toggle BL initialization flag

   Alfa r   Prescribe alpha
   CLI  r   Prescribe inviscid CL
   Cl   r   Prescribe CL
   ASeq rrr Prescribe a sequence of alphas
   CSeq rrr Prescribe a sequence of CLs

   SEQP     Toggle polar/Cp(x) sequence plot display
   CINC     Toggle  minimum Cp  inclusion in polar
   HINC     Toggle hinge moment inclusion in polar
   Pacc i   Toggle auto point accumulation to active polar
   PGET f   Read new polar from save file
   PWRT i   Write polar to save file
   PSUM     Show summary of stored polars
   PLIS i   List stored polar(s)
   PDEL i   Delete stored polar
   PSOR i   Sort stored polar
   PPlo ii. Plot stored polar(s)
   APlo ii. Plot stored airfoil(s) for each polar
   ASET i   Copy stored airfoil into current airfoil
   PREM ir. Remove point(s) from stored polar
   PNAM i   Change airfoil name of stored polar
   PPAX     Change polar plot axis limits

   RGET f   Read new reference polar from file
   RDEL i   Delete stored reference polar

   GRID     Toggle Cp vs x grid overlay
   CPIO     Toggle inviscid Cp overlay on viscous Cp(x)
   CREF     Toggle reference Cp data overlay
   FREF     Toggle reference CL,CD.. data display

   CPx      Plot Cp vs x
   CPV      Plot airfoil with pressure vectors (gee wiz)
  .VPlo     BL variable plots
  .ANNO     Annotate current plot
   HARD     Hardcopy current plot
   SIZE r   Change plot-object size
   CPMI r   Change minimum Cp axis annotation

   BL   i   Plot boundary layer velocity profiles
   BLC      Plot boundary layer velocity profiles at cursor
   BLWT r   Change velocity profile scale weight

   FMOM     Calculate flap hinge moment and forces
   FNEW rr  Set new flap hinge point
   VELS rr  Calculate velocity components at a point
   DUMP f   Output Ue,Dstar,Theta,Cf vs s,x,y to file
   CPWR f   Output x vs Cp to file
   CPMN     Report minimum surface Cp
   NAME s   Specify new airfoil name
   NINC     Increment name version number
)";

const char* const kVparMenu = R"(
   <cr>    Return to OPER menu
   SHOW    Display viscous parameters
   Xtr  rr Change trip positions Xtr/c
   N    r  Change critical amplification exponent Ncrit
   NT   r  Change Ncrit on Top
   NB   r  Change Ncrit on Bot
   Vacc r  Change Newton solution acceleration parameter
   Wake r  Change wake length/chord
   INIT    BL initialization flag toggle

   LAG     change lag equation constants
   GB      change G-beta constants A,B
   CTR     change initial transition-Ctau constants
   REST    restore BL calibration to baseline
)";

const char* const kGdesMenu = R"(
   <cr>     Return to Top Level
   !        Redo previous command

   GSET     Set buffer  airfoil <== current airfoil
   eXec     Set current airfoil <== buffer  airfoil
   SYMM     Toggle y-symmetry flag

   ADEG r   Rotate about origin (degrees)
   ARAD r   Rotate about origin (radians)
   Tran rr  Translate
   Scal r   Scale about origin
   LINS rr. Linearly-varying y scale
   DERO     Derotate (set chord line level)

   TGAP rr  Change trailing edge gap
   LERA rr  Change leading edge radius

   TCPL     Toggle thickness and camber plotting
   TFAC rr  Scale existing thickness and camber
   TSET rr  Set new thickness and camber
   HIGH rr  Move camber and thickness highpoints
  .CAMB     Modify camber shape directly or via loading

   BEND     Display structural properties of buffer airfoil

   Flap rrr Deflect trailing edge flap

   Modi     Modify contour via cursor
   SLOP     Toggle modified-contour slope matching flag

   CORN     Double point with cursor (set sharp corner)
   ADDP     Add    point with cursor or keyboard x,y
   MOVP     Move   point with cursor or keyboard x,y
   DELP     Delete point with cursor
   NMOV r   Move all points in surface-normal direction

   UNIT     Normalize buffer airfoil to unit chord
   Dist     Determine distance between 2 cursor points
   CLIS     List curvatures
   CPLO     Plot curvatures
   CANG     List panel corner angles
   CADD ri. Add points at corners exceeding angle threshold

   Plot     Replot buffer airfoil
   INPL     Replot buffer airfoil without scaling (in inches)
   Blow     Blowup plot region
   Rese     Reset plot scale and origin
   Wind     Plot window adjust via cursor and keys

   TSIZ r   Change tick-mark size
   TICK     Toggle node tick-mark plotting
   GRID     Toggle grid plotting
   GPAR     Toggle geometric parameter plotting
   Over f   Overlay disk file airfoil

   SIZE r   Change absolute plot-object size
  .ANNO     Annotate plot
   HARD     Hardcopy current plot

   NAME s   Specify new airfoil name
   NINC     Increment name version number
)";

const char* const kTypeTable =
    "\n Type   parameters held constant       varying      fixed   \n"
    R"( ----   ------------------------       -------   -----------
   1    M          , Re            ..   lift     chord, vel.
   2    M sqrt(CL) , Re sqrt(CL)   ..   vel.     chord, lift
   3    M          , Re CL         ..   chord    lift , vel.
)";

const char* const kPolarRule =
    " ==============================================================\n";

const char* const kVploMenu = R"(
   <cr>   Return to OPER menu

   H      Plot kinematic shape parameter
   DT     Plot top    side Dstar and Theta
   DB     Plot bottom side Dstar and Theta
   UE     Plot edge velocity
   CF     Plot skin friction coefficient
   CD     Plot dissipation coefficient
   N      Plot amplification ratio
   CT     Plot max shear coefficient
   RT     Plot Re_theta
   RTL    Plot log(Re_theta)

   DUMP f Write current plot variable to file
   OVER f Overlay current plot variable from file

   X  rrr Change x-axis limits
   Y  rrr Change y-axis limits on current plot
   XPLO   Toggle plotting  vs  x or s

   BLOW   Cursor blowup of current plot
   RESE   Reset to default x,y-axis limits
   SIZE r Change absolute plot-object size
  .ANNO   Annotate plot
   HARD   Hardcopy current plot

   GRID   Toggle grid plotting
   SYMB   Toggle node-symbol plotting
   LABE   Toggle label plotting
   CLIP   Toggle line-plot clipping
   FRPL   Toggle TS frequency plotting
)";

/// A real as gfortran's list-directed output writes it: 17 significant
/// digits, as F right-justified in 20 columns plus 5 blanks for magnitudes
/// from 0.1 to 1e16, as E right-justified in 25 columns otherwise.
std::string ListDirected(double value) {
    const double magnitude = std::fabs(value);
    char buffer[64];
    if (magnitude != 0.0 && (magnitude < 0.1 || magnitude >= 1e16)) {
        char mantissa[40];
        std::snprintf(mantissa, sizeof mantissa, "%.16e", value);
        std::string text = mantissa;
        const std::size_t e = text.find('e');
        const auto exponent =
            static_cast<int>(std::strtol(text.c_str() + e + 1, nullptr, 10));
        std::snprintf(buffer, sizeof buffer, "%sE%c%03d",
                      text.substr(0, e).c_str(), exponent < 0 ? '-' : '+',
                      std::abs(exponent));
        char field[64];
        std::snprintf(field, sizeof field, "%25s", buffer);
        return field;
    }
    /* Significant digits only: a leading "0." does not count */
    int integer_digits = magnitude == 0.0 ? 1 : 0;
    if (magnitude >= 1.0)
        integer_digits =
            static_cast<int>(std::floor(std::log10(magnitude))) + 1;
    const int decimals = std::max(0, 17 - integer_digits);
    std::snprintf(buffer, sizeof buffer, "%20.*f     ", decimals, value);
    return buffer;
}

/* Xplot11-style annotation menu, as XFOIL shows it */
const char* const kAnnoMenu =
    "\n   C haracters           |       \n"
    "   S lant characters     |       \n"
    "   M ath characters      |  plot \n"
    "   P oint symbol         |       \n"
    "   L ine                 |       \n"
    "   A rrow                |       \n"
    " \n"
    "   W idth of characters  | modify\n"
    "   T ype of point symbol |       \n"
    "\n";

const char* const kSymbolTable =
    "   0  square         7  Y         \n"
    "   1  circle         8  flipped Y \n"
    "   2  triangle       9  *         \n"
    "   3  +             10  flipped * \n"
    "   4  x             11  hourglass \n"
    "   5  diamond       12  bowtie    \n"
    "   6  yield sign    13  star      \n";

/// Number pairs from a text file, in blocks ended by a "999 999" line.
/// @return False if the file cannot be opened.
bool ReadPairBlocks(const std::string& path,
                    std::vector<std::vector<Point>>& blocks) {
    std::ifstream stream(path);
    if (!stream) return false;
    blocks.assign(1, {});
    std::string line;
    while (std::getline(stream, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        double a = 0.0;
        double b = 0.0;
        if (!(row >> a >> b)) continue;
        if (a == 999.0 && b == 999.0) {
            blocks.emplace_back();
            continue;
        }
        blocks.back().push_back({a, b});
    }
    if (blocks.back().empty() && blocks.size() > 1) blocks.pop_back();
    return true;
}

/// The page rectangle spanned by two points, widened to the page's aspect
/// ratio so a zoom keeps shapes undistorted.
Rect AspectRect(const Point& a, const Point& b) {
    Rect r{std::min(a.x, b.x), std::min(a.y, b.y), std::max(a.x, b.x),
           std::max(a.y, b.y)};
    const double aspect = Scene::kWidth / Scene::kHeight;
    const double cx = 0.5 * (r.x0 + r.x1);
    const double cy = 0.5 * (r.y0 + r.y1);
    double w = r.x1 - r.x0;
    double h = r.y1 - r.y0;
    if (w < aspect * h)
        w = aspect * h;
    else
        h = w / aspect;
    return {cx - 0.5 * w, cy - 0.5 * h, cx + 0.5 * w, cy + 0.5 * h};
}

void NotAvailable(std::ostream& out, const std::string& name) {
    out << ' ' << A(name, 4) << " command not available in Scirocco\n";
}

void NotRecognized(std::ostream& out, const std::string& name) {
    out << ' ' << A(name, 4)
        << " command not recognized.  Type a \"?\" for list\n";
}

/// Rotate points about the origin; positive angles turn the nose up, as
/// XFOIL's ROTATE does.
MatrixXd Rotated(const MatrixXd& x, double radians) {
    const double c = std::cos(radians);
    const double s = std::sin(radians);
    MatrixXd result(2, x.cols());
    result.row(0) = c * x.row(0) + s * x.row(1);
    result.row(1) = -s * x.row(0) + c * x.row(1);
    return result;
}

/// Upper and lower surface y at chordwise station @p xs.
std::pair<double, double> SurfaceY(const MatrixXd& x, double xs) {
    const int ile = LeadingEdgeIndex(x);
    auto along = [&](int begin, int end) {
        double best = 0.0;
        double best_distance = 1e300;
        for (int i = begin; i < end; ++i) {
            const double xa = x(0, i);
            const double xb = x(0, i + 1);
            if ((xs - xa) * (xs - xb) <= 0.0 && xa != xb) {
                const double w = (xs - xa) / (xb - xa);
                return (1.0 - w) * x(1, i) + w * x(1, i + 1);
            }
            const double d = std::fabs(xa - xs);
            if (d < best_distance) {
                best_distance = d;
                best = x(1, i);
            }
        }
        return best;
    };
    const double first = along(0, ile);
    const double second = along(ile, static_cast<int>(x.cols()) - 1);
    return {std::max(first, second), std::min(first, second)};
}

}  // namespace

Session::Session(std::istream& in, std::ostream& out, std::string version)
    : in_(in), out_(out), version_(std::move(version)) {}

Session::~Session() = default;

/* ------------------------------------------------------------------ */
/* Prompts                                                              */
/* ------------------------------------------------------------------ */

void Session::Prompt(const std::string& label, const std::string& kind) {
    out_ << '\n' << label << "   " << kind << ">  " << std::flush;
}

std::string Session::ReadLine() {
    std::string line;
    if (!std::getline(in_, line)) throw EndOfInput();
    if (!line.empty() && line.back() == '\r') line.pop_back();
    return line;
}

std::string Session::AskCommand(const std::string& label) {
    Prompt(label, "c");
    return ReadLine();
}

void Session::AskReal(const std::string& prompt, double& value) {
    for (;;) {
        Prompt(prompt, "r");
        const std::string line = Trim(ReadLine());
        if (line.empty()) return;
        const std::vector<double> values = Numbers(line);
        if (!values.empty()) {
            value = values[0];
            return;
        }
    }
}

void Session::AskInt(const std::string& prompt, int& value) {
    for (;;) {
        Prompt(prompt, "i");
        const std::string line = Trim(ReadLine());
        if (line.empty()) return;
        const std::vector<double> values = Numbers(line);
        if (!values.empty() && values[0] == std::floor(values[0])) {
            value = static_cast<int>(values[0]);
            return;
        }
    }
}

std::string Session::AskString(const std::string& prompt) {
    Prompt(prompt, "s");
    return Trim(ReadLine());
}

bool Session::AskYesNo(const std::string& prompt) {
    for (;;) {
        out_ << '\n' << prompt << " y/n>  " << std::flush;
        const std::string line = Trim(ReadLine());
        if (line.empty()) continue;
        const char c = static_cast<char>(std::toupper(line[0]));
        if (c == 'Y') return true;
        if (c == 'N') return false;
    }
}

bool Session::ConfirmOverwrite(const std::string& path, const char* refused) {
    if (!FileExists(path)) return true;
    out_ << "\n Output file exists.  Overwrite?  Y" << std::flush;
    const std::string line = Trim(ReadLine());
    out_ << '\n';
    if (!line.empty() && std::toupper(line[0]) == 'N') {
        if (refused) out_ << refused << '\n';
        return false;
    }
    return true;
}

/* ------------------------------------------------------------------ */
/* Top level                                                            */
/* ------------------------------------------------------------------ */

int Session::Run(const std::vector<std::string>& args) {
    out_ << '\n'
         << " ===================================================\n"
         << "  " << kProgram << " Version " << version_ << '\n'
         << "  Copyright (C) 2024-2026   Christopher A. Lupp\n"
         << '\n'
         << "  This software comes with ABSOLUTELY NO WARRANTY,\n"
         << "    subject to the Apache License, Version 2.0.\n"
         << '\n'
         << "  Caveat computor\n"
         << " ===================================================\n";
    try {
        if (!args.empty()) {
            out_ << '\n';
            Load(args[0]);
        }
        TopLevel();
    } catch (const EndOfInput&) {
        out_ << '\n';
    }
    out_ << std::flush;
    return 0;
}

void Session::TopLevel() {
    auto menu = [this] {
        out_ << kTopMenuHead << "   XYCM rr Change CM reference location, "
             << "currently " << F(xref_(0), 8, 5) << F(xref_(1), 8, 5)
             << "\n\n   BEND    Display structural properties of current "
                "airfoil\n\n"
             << "   PCOP    Set current-airfoil panel nodes directly from "
                "buffer airfoil points\n"
             << "   PANE    Set current-airfoil panel nodes (" << I(nodes_, 4)
             << " ) based on curvature" << kTopMenuTail;
    };
    menu();
    for (;;) {
        const Command command = Parse(AskCommand(" fauXFOIL"));
        if (command.name.empty()) continue;
        try {
            if (Is(command, {"QUIT", "Q"})) {
                return;
            } else if (command.name == "?") {
                menu();
            } else if (Is(command, {"OPER", "O"})) {
                Oper();
            } else if (command.name == "GDES") {
                Gdes();
            } else if (command.name == "NACA") {
                Naca(command.arg);
            } else if (command.name == "LOAD") {
                Load(command.arg);
            } else if (command.name == "SAVE") {
                Save(command.arg, true);
            } else if (command.name == "PSAV") {
                Save(command.arg, false);
            } else if (command.name == "REVE") {
                write_clockwise_ = !write_clockwise_;
                out_ << (write_clockwise_
                             ? " Airfoil will be written in clockwise order\n"
                             : " Airfoil will be written in counterclockwise "
                               "order\n");
            } else if (command.name == "NORM") {
                normalize_ = !normalize_;
                out_ << (normalize_ ? " Loaded airfoil will  be normalized\n"
                                    : " Loaded airfoil will not be "
                                      "normalized\n");
            } else if (command.name == "XYCM") {
                double x =
                    command.numbers.size() > 0 ? command.numbers[0] : xref_(0);
                double y =
                    command.numbers.size() > 1 ? command.numbers[1] : xref_(1);
                if (command.numbers.size() < 1)
                    AskReal("Enter new CM reference X", x);
                if (command.numbers.size() < 2)
                    AskReal("Enter new CM reference Y", y);
                xref_ = Vector2d(x, y);
            } else if (command.name == "PCOP") {
                if (buffer_.cols() == 0) {
                    out_ << " ***  No buffer airfoil available  ***\n";
                } else {
                    SetCurrentFromBuffer();
                    out_ << "\n Current airfoil nodes set from buffer "
                            "airfoil nodes ("
                         << I(buffer_.cols(), 4) << " )\n";
                }
            } else if (command.name == "PANE") {
                Pane();
            } else if (command.name == "PPAR") {
                Ppar();
            } else if (command.name == "PLOP") {
                Plop();
            } else if (command.name == "NAME") {
                Rename(command.arg, false);
            } else if (command.name == "NINC") {
                Rename("", true);
            } else if (Is(command, {"Z", "U"})) {
                ZoomCommand(command.name);
            } else if (Is(command, {"MDES", "QDES", "ISAV", "MSAV", "DELI",
                                    "INTE", "HALF", "BEND", "WDEF", "RDEF"})) {
                NotAvailable(out_, command.name);
            } else {
                NotRecognized(out_, command.name);
            }
        } catch (const EndOfInput&) {
            throw;
        } catch (const std::exception& error) {
            out_ << " *** " << error.what() << '\n';
        }
    }
}

void Session::Plop() {
    for (;;) {
        out_ << " ...............................................\n\n"
             << "  G raphics-enable flag:        " << (graphics_ ? 'T' : 'F')
             << '\n';
        const Command command =
            Parse(AskCommand("      Option, Value   (or <Return>) "));
        if (command.name.empty()) return;
        if (command.name[0] == 'G') graphics_ = !graphics_;
    }
}

/* ------------------------------------------------------------------ */
/* Geometry                                                             */
/* ------------------------------------------------------------------ */

void Session::PrintThicknessCamber(const MatrixXd& x) {
    const ThicknessCamber tc = GetThicknessCamber(x);
    out_ << " Max thickness = " << F(tc.thickness, 12, 6)
         << "  at x = " << F(tc.x_thickness, 7, 3) << '\n'
         << " Max camber    = " << F(tc.camber, 12, 6)
         << "  at x = " << F(tc.x_camber, 7, 3) << '\n';
}

void Session::PrintTrailingEdge() {
    const double gap = TrailingEdgeGap(foil_->GetFoilCoords());
    if (gap > 1e-9)
        out_ << "\n Blunt trailing edge.  Gap =" << F(gap, 9, 5) << '\n';
    else
        out_ << "\n Sharp trailing edge\n";
}

void Session::PrintPanelingUsed() {
    out_ << "\n Paneling parameters used...\n"
         << "   Number of panel nodes" << I(nodes_, 10) << '\n';
}

void Session::SetCurrentFromBuffer() {
    /* Panel nodes are the buffer points themselves (XFOIL's PCOP) */
    const int n = static_cast<int>(buffer_.cols());
    foil_ = std::make_unique<scirocco::Scirocco>(buffer_, n - 1);
    FoilState& state = foil_->GetState();
    state.foil.x = buffer_;
    state.foil.N = n;
    scirocco::RefreshPanelGeometry(state.foil);
    state.geom.name = name_;
    announce_vorticity_ = true;
    bl_initialized_ = false;
}

void Session::Pane() {
    if (buffer_.cols() == 0) {
        out_ << " ***  No buffer airfoil available  ***\n";
        return;
    }
    foil_ = std::make_unique<scirocco::Scirocco>(buffer_, nodes_ - 1);
    foil_->GetState().geom.name = name_;
    announce_vorticity_ = true;
    bl_initialized_ = false;
    PrintTrailingEdge();
    PrintPanelingUsed();
}

void Session::Naca(const std::string& arg) {
    std::string digits = Trim(arg);
    if (digits.empty()) {
        Prompt("Enter NACA 4 or 5-digit airfoil designation", "i");
        digits = Trim(ReadLine());
        if (digits.empty()) return;
    }
    digits = digits.substr(0, digits.find_first_of(" \t,"));
    if (digits.find_first_not_of("0123456789") != std::string::npos ||
        digits.size() > 5) {
        out_ << " *** Only 4 and 5-digit NACA airfoils are available\n";
        return;
    }
    if (digits.size() < 4) digits.insert(0, 4 - digits.size(), '0');
    scirocco::Geom geom;
    scirocco::NacaPoints(geom, digits);
    buffer_ = geom.xpoint;
    name_ = geom.name;
    PrintThicknessCamber(buffer_);
    out_ << "\n Buffer airfoil set using " << buffer_.cols() << " points\n";
    Pane();
}

void Session::Load(const std::string& arg) {
    std::string path = Trim(arg);
    if (path.empty()) path = AskString("Enter filename");
    if (path.empty()) return;
    AirfoilFile file;
    std::string error;
    if (!ReadAirfoilFile(path, file, error)) {
        if (error == "open")
            out_ << "\n File OPEN error.  Nonexistent file:  " << path
                 << "  \n";
        else
            out_ << "\n File READ error.  Unrecognized format:  " << path
                 << "  \n";
        out_ << " *** LOAD NOT COMPLETED ***\n";
        return;
    }
    std::string name = file.name;
    if (file.labeled) {
        out_ << "\n Labeled airfoil file.  Name:  " << A(name, 48) << '\n';
    } else {
        out_ << "\n Plain airfoil file\n";
        name = AskString("Enter airfoil name");
    }
    const bool ccw = IsCounterclockwise(file.x);
    MatrixXd points = ccw ? Reversed(file.x) : file.x;
    if (normalize_) points = Normalized(points);
    scirocco::ValidateCoordinates(points);

    buffer_ = points;
    name_ = name;
    out_ << "\n Number of input coordinate points:" << I(points.cols(), 4)
         << '\n'
         << (ccw ? " Counterclockwise ordering\n" : " Clockwise ordering\n");
    PrintThicknessCamber(buffer_);
    const Vector2d le = LeadingEdge(buffer_);
    const Vector2d te =
        0.5 * (buffer_.col(0) + buffer_.col(buffer_.cols() - 1));
    out_ << "\n  LE  x,y  =" << F(le(0), 10, 5) << F(le(1), 10, 5)
         << "  |   Chord =" << F(Chord(buffer_), 10, 5) << '\n'
         << "  TE  x,y  =" << F(te(0), 10, 5) << F(te(1), 10, 5) << "  |\n";
    SetCurrentFromBuffer();
    out_ << "\n Current airfoil nodes set from buffer airfoil nodes ("
         << I(buffer_.cols(), 4) << " )\n";
}

void Session::Save(const std::string& arg, bool labeled) {
    if (!foil_) {
        out_ << " ***  No current airfoil to save  ***\n";
        return;
    }
    std::string path = Trim(arg);
    if (path.empty()) path = AskString("Enter output filename");
    if (path.empty()) return;
    if (!ConfirmOverwrite(path, " Current airfoil not saved.")) return;
    const MatrixXd& x = foil_->GetFoilCoords();
    const MatrixXd points = write_clockwise_ ? x : Reversed(x);
    if (!WriteAirfoilFile(path, labeled ? &name_ : nullptr, points))
        out_ << " *** File write error: " << path << '\n';
}

void Session::Rename(const std::string& arg, bool increment) {
    if (!increment && !Trim(arg).empty()) {
        name_ = Trim(arg);
    } else {
        std::string fallback = name_;
        if (increment) {
            /* "NAME" -> "NAME [1]" -> "NAME [2]" ... */
            const std::size_t open = fallback.rfind(" [");
            int version = 0;
            if (open != std::string::npos && fallback.back() == ']') {
                version = static_cast<int>(
                    std::strtol(fallback.c_str() + open + 2, nullptr, 10));
                fallback = fallback.substr(0, open);
            }
            fallback += " [" + std::to_string(version + 1) + "]";
        }
        out_ << "\n Enter airfoil name or <return> for default:  "
             << A(fallback, 48) << '\n'
             << std::flush;
        const std::string line = Trim(ReadLine());
        name_ = line.empty() ? fallback : line;
    }
    if (foil_) foil_->GetState().geom.name = name_;
}

bool Session::BufferMatchesCurrent() const {
    if (!foil_) return buffer_.cols() == 0;
    const MatrixXd& x = foil_->GetFoilCoords();
    return x.cols() == buffer_.cols() && x.isApprox(buffer_, 1e-12);
}

void Session::Ppar() {
    auto show = [this] {
        out_ << "\n    Present paneling parameters...\n"
             << "  N  i   Number of panel nodes" << I(nodes_, 10) << '\n';
    };
    show();
    bool changed = false;
    for (;;) {
        const Command command =
            Parse(AskCommand("Change what ? (<cr> if nothing else)"));
        if (command.name.empty()) {
            if (!changed) return;
            changed = false;
            if (buffer_.cols() == 0) return;
            foil_ = std::make_unique<scirocco::Scirocco>(buffer_, nodes_ - 1);
            foil_->GetState().geom.name = name_;
            announce_vorticity_ = true;
            bl_initialized_ = false;
            PrintTrailingEdge();
            const MatrixXd& x = foil_->GetFoilCoords();
            int node = 0;
            const double angle = MaxCornerAngle(x, &node);
            out_ << "\n Maximum panel corner angle =" << F(angle, 7, 3)
                 << "   at  i,x,y  =" << I(x.cols() - node, 4)
                 << F(x(0, node), 9, 4) << F(x(1, node), 9, 4) << '\n';
            show();
        } else if (command.name == "N") {
            int n = nodes_;
            if (!command.numbers.empty())
                n = static_cast<int>(command.numbers[0]);
            else
                AskInt("Enter number of panel nodes", n);
            if (n >= 3 && n - 1 <= scirocco::kMaxPanels) {
                changed = changed || n != nodes_;
                nodes_ = n;
            }
        } else if (command.name == "?") {
            show();
        } else {
            out_ << "\n ***  Input not recognized  ***\n";
            show();
        }
    }
}

void Session::Gdes() {
    if (buffer_.cols() == 0) {
        out_ << " ***  No buffer airfoil available  ***\n";
        return;
    }
    out_ << "\n You are working with the buffer airfoil\n";
    PlotAirfoil();
    Command previous;
    for (;;) {
        Command command = Parse(AskCommand(".GDES"));
        if (command.name == "!") command = previous;
        if (command.name.empty()) break;
        previous = command;
        const auto& v = command.numbers;
        if (command.name == "?") {
            out_ << kGdesMenu;
        } else if (command.name == "GSET") {
            if (foil_) buffer_ = foil_->GetFoilCoords();
            PlotAirfoil();
        } else if (Is(command, {"EXEC", "X"})) {
            SetCurrentFromBuffer();
            out_ << "\n Current airfoil nodes set from buffer airfoil nodes ("
                 << I(buffer_.cols(), 4) << " )\n";
        } else if (Is(command, {"ADEG", "ARAD"})) {
            double angle = v.empty() ? 0.0 : v[0];
            if (v.empty())
                AskReal(command.name == "ADEG" ? "Enter angle change (deg)"
                                               : "Enter angle change (rad)",
                        angle);
            if (command.name == "ADEG") angle *= kPi / 180.0;
            buffer_ = Rotated(buffer_, angle);
            PrintThicknessCamber(buffer_);
            PlotAirfoil();
        } else if (Is(command, {"TRAN", "T"})) {
            double dx = v.size() > 0 ? v[0] : 0.0;
            double dy = v.size() > 1 ? v[1] : 0.0;
            if (v.size() < 1) AskReal("Enter delta(x)", dx);
            if (v.size() < 2) AskReal("Enter delta(y)", dy);
            buffer_.row(0).array() += dx;
            buffer_.row(1).array() += dy;
            PlotAirfoil();
        } else if (Is(command, {"SCAL", "S"})) {
            double scale = v.empty() ? 1.0 : v[0];
            if (v.empty()) AskReal("Enter scale factor", scale);
            buffer_ *= scale;
            PrintThicknessCamber(buffer_);
            PlotAirfoil();
        } else if (command.name == "UNIT") {
            buffer_ = Normalized(buffer_);
            PrintThicknessCamber(buffer_);
            PlotAirfoil();
        } else if (command.name == "DERO") {
            auto chord_angle = [](const MatrixXd& x) {
                const Vector2d le = LeadingEdge(x);
                const Vector2d te = 0.5 * (x.col(0) + x.col(x.cols() - 1));
                return std::atan2(te(1) - le(1), te(0) - le(0)) * 180.0 / kPi;
            };
            const double before = chord_angle(buffer_);
            FoilState scratch;
            scratch.geom.xpoint = buffer_;
            scratch.geom.npoint = static_cast<int>(buffer_.cols());
            scirocco::GeomDerotate(scratch);
            buffer_ = scratch.geom.xpoint;
            const double angle = before - chord_angle(buffer_);
            out_ << "\nRotating buffer airfoil by " << F(angle, 8, 3)
                 << " deg.\n";
            PrintThicknessCamber(buffer_);
            PlotAirfoil();
        } else if (Is(command, {"FLAP", "F"})) {
            double x = v.size() > 0 ? v[0] : 0.0;
            double y = v.size() > 1 ? v[1] : 0.0;
            double eta = v.size() > 2 ? v[2] : 0.0;
            if (v.size() < 1) AskReal("Enter flap hinge x location", x);
            const auto [yu, yl] = SurfaceY(buffer_, x);
            if (v.size() < 2) {
                out_ << "\n  Top    surface:  y =" << F(yu, 8, 4)
                     << "     y/t = 1.0\n"
                     << "  Bottom surface:  y =" << F(yl, 8, 4)
                     << "     y/t = 0.0\n";
                AskReal("Enter flap hinge y location (or 999 to specify y/t)",
                        y);
            }
            if (y == 999.0) {
                double fraction = 0.5;
                AskReal("Enter flap hinge relative y/t", fraction);
                y = yl + fraction * (yu - yl);
            }
            out_ << "\n Flap hinge: x,y =" << F(x, 9, 5) << F(y, 9, 5) << '\n';
            if (v.size() < 3)
                AskReal("Enter flap deflection in degrees (+ down)", eta);
            FoilState scratch;
            scratch.geom.xpoint = buffer_;
            scratch.geom.npoint = static_cast<int>(buffer_.cols());
            scirocco::GeomFlap(scratch, Vector2d(x, y), eta);
            buffer_ = scratch.geom.xpoint;
            out_ << '\n';
            PrintThicknessCamber(buffer_);
            PlotAirfoil();
        } else if (command.name == "NAME") {
            Rename(command.arg, false);
        } else if (command.name == "NINC") {
            Rename("", true);
        } else if (Is(command, {"PLOT", "P"})) {
            airfoil_options_.inches = false;
            PlotAirfoil();
        } else if (command.name == "GRID") {
            airfoil_options_.grid = !airfoil_options_.grid;
            PlotAirfoil();
        } else if (command.name == "TICK") {
            airfoil_options_.ticks = !airfoil_options_.ticks;
            PlotAirfoil();
        } else if (command.name == "GPAR") {
            airfoil_options_.parameters = !airfoil_options_.parameters;
            PlotAirfoil();
        } else if (command.name == "HARD") {
            Hardcopy();
        } else if (Is(command, {"BLOW", "B"})) {
            GdesBlowup();
        } else if (Is(command, {"RESE", "R"})) {
            airfoil_options_.window_set = false;
            airfoil_options_.inches = false;
            PlotAirfoil();
        } else if (Is(command, {"WIND", "W"})) {
            GdesWindow();
        } else if (command.name == "INPL") {
            airfoil_options_.inches = true;
            airfoil_options_.window_set = false;
            PlotAirfoil();
        } else if (command.name == "TSIZ") {
            TickSize(v);
        } else if (Is(command, {"OVER", "O"})) {
            GdesOverlay(command.arg);
        } else if (command.name == "CPLO") {
            Curvatures();
        } else if (Is(command, {"DIST", "D"})) {
            Distance();
        } else if (command.name == "SIZE") {
            SizeCommand(gdes_size_, v);
        } else if (command.name == "ANNO") {
            Annotate();
        } else if (Is(command, {"Z", "U"})) {
            ZoomCommand(command.name);
        } else if (Is(command, {"SYMM", "LINS", "TGAP", "LERA", "TCPL", "TFAC",
                                "TSET", "HIGH", "CAMB", "BEND", "MODI", "M",
                                "SLOP", "CORN", "ADDP", "MOVP", "DELP", "NMOV",
                                "CLIS", "CANG", "CADD"})) {
            NotAvailable(out_, command.name);
        } else {
            out_ << " Command " << A(command.name, 4)
                 << " not recognized.  Type a \" ? \" for list.\n";
        }
    }
    if (!BufferMatchesCurrent())
        out_ << "\n Buffer airfoil is not identical to current airfoil\n";
}

/* ------------------------------------------------------------------ */
/* OPER                                                                 */
/* ------------------------------------------------------------------ */

std::string Session::OperLabel() const {
    std::string label = viscous_ ? ".OPERv" : ".OPERi";
    if (pacc_) label += "a";
    return label;
}

bool Session::PolarLocked() {
    if (!pacc_) return false;
    out_ << "\n * Polar is being accumulated.\n"
         << " * Cannot change its parameters in midstream.\n";
    return true;
}

void Session::PrintMachRe() const {
    out_ << "\n M  =" << F(mach_, 10, 4) << '\n'
         << " Re =" << E(re_, 12, 4) << '\n';
}

void Session::PrintSonic(bool extra_blank) const {
    if (mach_ <= 0.0) return;
    constexpr double gamma = 1.4;
    const double m2 = mach_ * mach_;
    const double ratio =
        (1.0 + 0.5 * (gamma - 1.0) * m2) / (1.0 + 0.5 * (gamma - 1.0));
    const double cp_star =
        2.0 / (gamma * m2) * (std::pow(ratio, gamma / (gamma - 1.0)) - 1.0);
    const double q_star = std::sqrt(ratio) / mach_;
    out_ << "\n Sonic Cp =" << F(cp_star, 10, 2)
         << "      Sonic Q/Qinf =" << F(q_star, 10, 3) << '\n';
    if (extra_blank) out_ << '\n';
}

void Session::Oper() {
    for (;;) {
        Command command = Parse(AskCommand(OperLabel()));
        if (command.name.empty()) return;
        if (command.name == "!") {
            if (last_solve_.empty()) continue;
            command = Parse(last_solve_);
        }
        const auto& v = command.numbers;
        try {
            if (command.name == "?") {
                out_ << kOperMenu;
            } else if (Is(command, {"VISC", "V"})) {
                if (PolarLocked()) continue;
                viscous_ = !viscous_;
                if (viscous_) {
                    if (!v.empty())
                        re_ = v[0];
                    else if (re_ <= 0.0)
                        AskReal("Enter Reynolds number", re_);
                    PrintMachRe();
                }
            } else if (command.name == "RE") {
                if (viscous_ && PolarLocked()) continue;
                if (!v.empty()) {
                    re_ = v[0];
                } else {
                    out_ << "\n Currently...\n\n Re =" << G(re_, 12, 4) << '\n';
                    AskReal("Enter new Reynolds number", re_);
                }
            } else if (Is(command, {"MACH", "M"})) {
                if (PolarLocked()) continue;
                double mach = mach_;
                if (!v.empty()) {
                    mach = v[0];
                } else {
                    out_ << "\n Currently...\n\n M  =" << F(mach_, 10, 4)
                         << '\n';
                    AskReal("Enter Mach number", mach);
                }
                if (mach < 0.0 || mach >= 1.0) {
                    out_ << " *** Mach number must satisfy 0 <= M < 1\n";
                    continue;
                }
                mach_ = mach;
                PrintSonic(viscous_ && !v.empty());
            } else if (Is(command, {"TYPE", "T"})) {
                if (PolarLocked()) continue;
                int type = 1;
                if (!v.empty()) {
                    type = static_cast<int>(v[0]);
                } else {
                    out_ << kTypeTable;
                    AskInt("Enter type of Mach,Re variation with CL", type);
                }
                if (type != 1)
                    out_ << " *** Only type 1 (fixed M and Re) is available "
                            "in Scirocco\n";
                PrintMachRe();
            } else if (command.name == "ITER") {
                if (!v.empty()) {
                    bl_.niglob = std::max(1, static_cast<int>(v[0]));
                } else {
                    out_ << " Current iteration limit:" << I(bl_.niglob, 12)
                         << '\n';
                    AskInt("Enter new iteration limit", bl_.niglob);
                    bl_.niglob = std::max(1, bl_.niglob);
                }
            } else if (command.name == "INIT") {
                bl_initialized_ = !bl_initialized_;
                out_ << (bl_initialized_
                             ? " BLs are assumed to be initialized\n"
                             : " BLs will be initialized on next point\n");
            } else if (Is(command, {"ALFA", "A"})) {
                double alpha = foil_ ? foil_->GetAlpha() : 0.0;
                if (!v.empty())
                    alpha = v[0];
                else
                    AskReal("Enter angle of attack (deg)", alpha);
                last_solve_ = "ALFA " + std::to_string(alpha);
                SolvePoint(Target::Alpha, alpha);
            } else if (Is(command, {"CL", "C"})) {
                double cl = foil_ ? foil_->GetCl() : 0.0;
                if (!v.empty())
                    cl = v[0];
                else
                    AskReal("Enter lift coefficient", cl);
                last_solve_ = "CL " + std::to_string(cl);
                SolvePoint(Target::Cl, cl);
            } else if (command.name == "CLI") {
                double cl = foil_ ? foil_->GetCl() : 0.0;
                if (!v.empty())
                    cl = v[0];
                else
                    AskReal("Enter inviscid lift coefficient", cl);
                last_solve_ = "CLI " + std::to_string(cl);
                SolvePoint(Target::InviscidCl, cl);
            } else if (Is(command, {"ASEQ", "AS", "CSEQ", "CS"})) {
                const bool alpha = command.name[0] == 'A';
                std::vector<double> range = v;
                range.resize(3, 0.0);
                const char* const prompts[2][3] = {
                    {"Enter first CL value", "Enter last  CL value",
                     "Enter CL increment  "},
                    {"Enter first alfa value (deg)",
                     "Enter last  alfa value (deg)",
                     "Enter alfa increment   (deg)"}};
                for (std::size_t k = v.size(); k < 3; ++k)
                    AskReal(prompts[alpha ? 1 : 0][k], range[k]);
                std::ostringstream redo;
                redo.precision(17);
                redo << (alpha ? "ASEQ " : "CSEQ ") << range[0] << ' '
                     << range[1] << ' ' << range[2];
                last_solve_ = redo.str();
                Sequence(alpha ? Target::Alpha : Target::Cl, range);
            } else if (Is(command, {"PACC", "P"})) {
                TogglePacc();
            } else if (command.name == "PGET") {
                ReadPolar(command.arg);
            } else if (command.name == "PWRT") {
                WritePolar(v);
            } else if (command.name == "PSUM") {
                PolarSummary();
            } else if (command.name == "PLIS") {
                ListPolars(v.empty() ? 0 : static_cast<int>(v[0]));
            } else if (command.name == "PDEL") {
                DeletePolar(v);
            } else if (command.name == "CPWR") {
                WriteCp(command.arg);
            } else if (command.name == "DUMP") {
                WriteDump(command.arg);
            } else if (command.name == "CPMN") {
                MinimumCp();
            } else if (command.name == "NAME") {
                Rename(command.arg, false);
            } else if (command.name == "NINC") {
                Rename("", true);
            } else if (command.name == "VPAR") {
                Vpar();
            } else if (Is(command, {"CPX", "CP"})) {
                if (HaveSolution()) PlotCp();
            } else if (command.name == "CPV") {
                if (HaveSolution()) PlotCpVectors();
            } else if (Is(command, {"VPLO", "VP"})) {
                Vplo();
            } else if (Is(command, {"PPLO", "PP"})) {
                if (!polars_.empty()) PlotPolars();
            } else if (command.name == "HARD") {
                Hardcopy();
            } else if (command.name == "GRID") {
                cp_options_.grid = !cp_options_.grid;
                out_ << (cp_options_.grid ? " Cp grid overlay enabled\n"
                                          : " Cp grid overlay disabled\n");
                if (plot_kind_ == PlotKind::Cp ||
                    plot_kind_ == PlotKind::CpSequence)
                    Replot();
            } else if (command.name == "CPIO") {
                cp_options_.inviscid = !cp_options_.inviscid;
                out_ << (cp_options_.inviscid
                             ? " Inviscid Cp overlay enabled\n"
                             : " Inviscid Cp overlay disabled\n");
                if (plot_kind_ == PlotKind::Cp) Replot();
            } else if (command.name == "SEQP") {
                sequence_polar_ = !sequence_polar_;
                out_ << (sequence_polar_
                             ? " Polar will be plotted during point sequence\n"
                             : " Cp(x) will be plotted during point "
                               "sequence\n");
            } else if (command.name == "CPMI") {
                double cp_min = cp_options_.cp_min;
                if (!v.empty()) {
                    cp_min = v[0];
                } else {
                    out_ << " Current CPmin = " << ListDirected(cp_min) << '\n';
                    AskReal("Enter new CPmin", cp_min);
                }
                cp_options_.cp_min = cp_min;
                if (plot_kind_ == PlotKind::Cp ||
                    plot_kind_ == PlotKind::CpSequence)
                    Replot();
            } else if (Is(command, {"Z", "U"})) {
                ZoomCommand(command.name);
            } else if (command.name == "SIZE") {
                SizeCommand(oper_size_, v);
            } else if (command.name == "ANNO") {
                Annotate();
            } else if (command.name == "BL") {
                Profiles(v);
            } else if (command.name == "BLC") {
                ProfilesAtCursor();
            } else if (command.name == "BLWT") {
                ProfileWeight(v);
            } else if (command.name == "CREF") {
                cp_reference_on_ = !cp_reference_on_;
                out_ << (cp_reference_on_
                             ? " Reference Cp plotting enabled\n"
                             : " Reference Cp plotting disabled\n");
            } else if (command.name == "FREF") {
                force_reference_on_ = !force_reference_on_;
                out_ << (force_reference_on_
                             ? " Reference force plotting enabled\n"
                             : " Reference force plotting disabled\n");
            } else if (command.name == "RGET") {
                ReadReferencePolar(command.arg);
            } else if (command.name == "RDEL") {
                DeleteReferencePolar(v);
            } else if (command.name == "PPAX") {
                PolarAxesCommand();
            } else if (Is(command, {"APLO", "AP"})) {
                if (!polars_.empty()) {
                    std::vector<int> which;
                    which.reserve(v.size());
                    for (double k : v) which.push_back(static_cast<int>(k));
                    Show(PolarAirfoilsPlot(polars_, which),
                         PlotKind::PolarAirfoils);
                }
            } else if (Is(command, {"CINC", "HINC", "PSOR", "ASET", "PREM",
                                    "PNAM", "FMOM", "FNEW", "VELS"})) {
                NotAvailable(out_, command.name);
            } else {
                NotRecognized(out_, command.name);
            }
        } catch (const EndOfInput&) {
            throw;
        } catch (const std::exception& error) {
            out_ << " *** " << error.what() << '\n';
        }
    }
}

void Session::ShowVpar() const {
    out_ << "\n Xtr/c     =" << F(1.0, 8, 4) << "    top    side\n"
         << " Xtr/c     =" << F(1.0, 8, 4) << "    bottom side\n"
         << " NcritT    =" << F(bl_.ncrit, 8, 2) << "   ("
         << F(TurbulenceLevel(bl_.ncrit), 6, 3) << " % turb. level )\n"
         << " NcritB    =" << F(bl_.ncrit, 8, 2) << "   ("
         << F(TurbulenceLevel(bl_.ncrit), 6, 3) << " % turb. level )\n"
         << " Vacc      =" << F(vacc_, 8, 4) << '\n'
         << " WakeL/c   =" << F(wake_length_, 8, 3) << '\n'
         << '\n'
         << " Klag  =" << F(bl_.SlagK, 8, 4)
         << "     Uxwt  =" << F(bl_.Cuq, 8, 2)
         << "       Kdl =" << F(bl_.Dlr, 8, 4) << '\n'
         << " A     =" << F(bl_.GA, 8, 4) << "     B     =" << F(bl_.GB, 8, 4)
         << "       KCt =" << F(0.5 / (bl_.GA * bl_.GA * bl_.GB), 8, 5) << '\n'
         << " CtiniK=" << F(bl_.CtauC, 8, 4)
         << "     CtiniX=" << F(bl_.CtauE, 8, 4) << '\n';
}

void Session::Vpar() {
    ShowVpar();
    for (;;) {
        const Command command = Parse(AskCommand("..VPAR"));
        if (command.name.empty()) return;
        const auto& v = command.numbers;
        if (command.name == "?") {
            out_ << kVparMenu;
        } else if (command.name == "SHOW") {
            ShowVpar();
        } else if (Is(command, {"XTR", "X"})) {
            if (PolarLocked()) continue;
            double top = v.size() > 0 ? v[0] : 1.0;
            double bottom = v.size() > 1 ? v[1] : 1.0;
            if (v.size() < 1) AskReal("Enter top    side Xtrip/c", top);
            if (v.size() < 2) AskReal("Enter bottom side Xtrip/c", bottom);
            if (top < 1.0 || bottom < 1.0)
                out_ << " *** Forced transition is not available in "
                        "Scirocco; Xtr/c stays 1.0\n";
        } else if (Is(command, {"N", "NT", "NB"})) {
            if (PolarLocked()) continue;
            double ncrit = bl_.ncrit;
            if (!v.empty())
                ncrit = v[0];
            else
                AskReal("Enter critical amplification ratio", ncrit);
            bl_.ncrit = ncrit;
            if (command.name != "N")
                out_ << " *** Scirocco uses one Ncrit for both sides; "
                        "NcritT and NcritB both set\n";
        } else if (command.name == "VACC") {
            if (!v.empty())
                vacc_ = v[0];
            else
                AskReal("Enter Newton acceleration parameter", vacc_);
            out_ << " *** Vacc has no effect in Scirocco\n";
        } else if (command.name == "WAKE") {
            if (!v.empty())
                wake_length_ = v[0];
            else
                AskReal("Enter wake length/chord", wake_length_);
        } else if (command.name == "INIT") {
            bl_initialized_ = !bl_initialized_;
            out_ << (bl_initialized_
                         ? " BLs are assumed to be initialized\n"
                         : " BLs will be initialized on next point\n");
        } else if (command.name == "LAG") {
            AskReal("Enter shear lag constant Klag", bl_.SlagK);
            AskReal("Enter shear lag UxEQ weight Uxwt", bl_.Cuq);
        } else if (command.name == "GB") {
            AskReal("Enter G-beta constant A", bl_.GA);
            AskReal("Enter G-beta constant B", bl_.GB);
        } else if (command.name == "CTR") {
            AskReal("Enter initial-Ctau constant CtiniK", bl_.CtauC);
            AskReal("Enter initial-Ctau exponent CtiniX", bl_.CtauE);
        } else if (command.name == "REST") {
            const scirocco::Param defaults;
            bl_.SlagK = defaults.SlagK;
            bl_.Cuq = defaults.Cuq;
            bl_.Dlr = defaults.Dlr;
            bl_.GA = defaults.GA;
            bl_.GB = defaults.GB;
            bl_.CtauC = defaults.CtauC;
            bl_.CtauE = defaults.CtauE;
        } else {
            NotRecognized(out_, command.name);
        }
    }
}

/* ------------------------------------------------------------------ */
/* Analysis                                                             */
/* ------------------------------------------------------------------ */

void Session::ApplyParameters(FoilState& state) const {
    state.param.niglob = bl_.niglob;
    state.param.ncrit = bl_.ncrit;
    state.param.SlagK = bl_.SlagK;
    state.param.Cuq = bl_.Cuq;
    state.param.Dlr = bl_.Dlr;
    state.param.GA = bl_.GA;
    state.param.GB = bl_.GB;
    state.param.CtauC = bl_.CtauC;
    state.param.CtauE = bl_.CtauE;
    state.geom.wakelen = wake_length_;
    state.geom.xref = xref_;
    state.geom.name = name_;
}

void Session::OnIteration(const FoilState& state, int iteration) {
    /* Iteration 0 is the initial guess; report each Newton step after it */
    if (iteration == 0) return;
    ++iteration_;
    const Transition top = FindTransition(state, 1);
    const Transition bottom = FindTransition(state, 0);
    out_ << '\n';
    out_ << " Side 1 " << (top.free ? " free " : "forced")
         << " transition at x/c = " << F(top.x, 7, 4) << I(top.station, 5)
         << '\n';
    out_ << " Side 2 " << (bottom.free ? " free " : "forced")
         << " transition at x/c = " << F(bottom.x, 7, 4) << I(bottom.station, 5)
         << '\n';
    out_ << '\n';

    /* XFOIL's convergence measure: Newton changes normalized by the
       variables (amplification by 10, edge velocity by 0.25 Vinf) */
    const auto& U = state.glob.U;
    const auto& dU = state.glob.dU;
    double sum = 0.0;
    int count = 0;
    double max_change = 0.0;
    char max_variable = 'C';
    int max_station = 0;
    int max_side = 1;
    if (dU.cols() == U.cols()) {
        for (int si = 0; si < 3; ++si) {
            const auto& stations = state.vsol.Is[si];
            for (std::size_t k = 0; k < stations.size(); ++k) {
                const int j = stations[k];
                const bool turbulent =
                    static_cast<std::size_t>(j) < state.vsol.turb.size() &&
                    state.vsol.turb[j];
                auto ratio = [](double d, double value) {
                    return value != 0.0 ? d / value : 0.0;
                };
                const double changes[4] = {
                    turbulent ? ratio(dU(2, j), U(2, j)) : dU(2, j) / 10.0,
                    ratio(dU(0, j), U(0, j)), ratio(dU(1, j), U(1, j)),
                    dU(3, j) / (0.25 * state.oper.Vinf)};
                const char letters[4] = {'C', 'T', 'D', 'U'};
                for (int m = 0; m < 4; ++m) {
                    sum += changes[m] * changes[m];
                    ++count;
                    if (std::fabs(changes[m]) > std::fabs(max_change)) {
                        max_change = changes[m];
                        max_variable = letters[m];
                        max_side = si == 1 ? 1 : 2;
                        max_station =
                            static_cast<int>(k) + 2 +
                            (si == 2 ? static_cast<int>(state.vsol.Is[0].size())
                                     : 0);
                    }
                }
            }
        }
    }
    const double rms = count > 0 ? std::sqrt(sum / count) : 0.0;
    out_ << I(iteration_, 4) << "   rms: " << E(rms, 10, 4)
         << "   max: " << E(max_change, 10, 4) << "   " << max_variable
         << " at " << I(max_station, 4) << I(max_side, 3);
    if (state.glob.relax < 1.0) out_ << "   RLX:" << F(state.glob.relax, 6, 3);
    out_ << '\n';
    out_ << "       a =" << F(state.oper.alpha, 7, 3)
         << "      CL =" << F(state.post.cl, 8, 4) << '\n';
    out_ << "      Cm =" << F(state.post.cm, 8, 4)
         << "     CD =" << F(state.post.cd, 9, 5)
         << "   =>   CDf =" << F(state.post.cdf, 9, 5)
         << "    CDp =" << F(state.post.cd - state.post.cdf, 9, 5) << '\n'
         << std::flush;
}

void Session::SolvePoint(Target target, double value) {
    if (!foil_) {
        out_ << " ***  No airfoil available  ***\n";
        return;
    }
    if (announce_vorticity_) {
        out_ << " Calculating unit vorticity distributions ...\n";
        announce_vorticity_ = false;
    }
    FoilState& state = foil_->GetState();
    ApplyParameters(state);

    double alpha = state.oper.alpha;
    if (target == Target::Alpha) alpha = value;
    if (target == Target::InviscidCl && viscous_) {
        /* Angle at which the inviscid solution gives this CL */
        foil_->SetOper(alpha, -1.0, 0.0, mach_);
        state.oper.givencl = true;
        state.oper.cltgt = value;
        if (!foil_->Solve()) {
            out_ << " *** Inviscid CL not reached\n";
            return;
        }
        alpha = foil_->GetAlpha();
        target = Target::Alpha;
    }
    const bool given_cl = target != Target::Alpha;

    bool converged = false;
    if (viscous_) {
        out_ << " Calculating wake trajectory ...\n"
             << " Calculating source influence matrix ...\n\n"
             << " Solving BL system ...\n";
        const bool cold = !bl_initialized_;
        if (cold)
            out_ << "\n Initializing BL ...\n"
                 << "    side            1  ...\n"
                 << "    side            2  ...\n";
        foil_->SetOper(alpha, -1.0, re_, mach_);
        state.oper.givencl = given_cl;
        if (given_cl) state.oper.cltgt = value;
        state.param.warm_start = !cold;
        iteration_ = 0;
        state.glob.iteration_observer = [this](const FoilState& s, int it) {
            OnIteration(s, it);
        };
        try {
            converged = foil_->Solve();
        } catch (const std::runtime_error&) {
            converged = false;
        }
        FoilState& after = foil_->GetState();
        after.glob.iteration_observer = nullptr;
        after.param.warm_start = true;
        if (converged)
            bl_initialized_ = true;
        else
            out_ << "\n VISCAL:  Convergence failed\n";
    } else {
        foil_->SetOper(alpha, -1.0, 0.0, mach_);
        state.oper.givencl = given_cl;
        if (given_cl) state.oper.cltgt = value;
        converged = foil_->Solve();
        const FoilState& after = foil_->GetState();
        /* XFOIL shows these only on the Cp plot; print them when no plot
           window shows it */
        if (display_ == nullptr || !graphics_)
            out_ << "       a =" << F(after.oper.alpha, 7, 3)
                 << "      CL =" << F(after.post.cl, 8, 4) << '\n'
                 << "      Cm =" << F(after.post.cm, 8, 4)
                 << "    CDp =" << F(after.post.cdpi, 9, 5) << '\n';
        if (!converged) out_ << " *** CL target not reached\n";
    }
    if (converged && pacc_) AddPolarPoint();
    out_ << std::flush;
    PlotPoint();
}

void Session::Sequence(Target target, const std::vector<double>& range) {
    const double first = range[0];
    const double last = range[1];
    double step = range[2];
    if ((last - first) * step < 0.0) step = -step;
    const int count =
        step == 0.0
            ? 1
            : static_cast<int>(std::lround(std::fabs((last - first) / step))) +
                  1;
    in_sequence_ = true;
    sequence_.clear();
    try {
        for (int k = 0; k < count; ++k) SolvePoint(target, first + k * step);
    } catch (...) {
        in_sequence_ = false;
        throw;
    }
    in_sequence_ = false;
}

/* ------------------------------------------------------------------ */
/* Polars                                                               */
/* ------------------------------------------------------------------ */

void Session::WritePolarFile(const Polar& polar) {
    std::ofstream stream(polar.save_file);
    stream << FormatPolar(polar, kProgram, version_);
}

void Session::TogglePacc() {
    if (pacc_) {
        pacc_ = false;
        out_ << "\n Polar accumulation disabled\n";
        return;
    }
    Polar polar;
    polar.name = name_;
    polar.re = viscous_ ? re_ : 0.0;
    polar.mach = mach_;
    polar.ncrit_top = polar.ncrit_bot = bl_.ncrit;
    if (foil_) {
        const MatrixXd& x = foil_->GetFoilCoords();
        for (Eigen::Index i = 0; i < x.cols(); ++i) {
            polar.airfoil_x.push_back(x(0, i));
            polar.airfoil_y.push_back(x(1, i));
        }
    }
    const int index = static_cast<int>(polars_.size()) + 1;
    out_ << "\n Polar" << I(index, 3) << " newly created for accumulation\n"
         << " Airfoil archived with polar: " << A(name_, 48) << '\n';

    const std::string save =
        AskString("Enter  polar save filename  OR  <return> for no file");
    if (save.empty()) {
        out_ << "\n Polar save file will NOT be written\n";
    } else if (FileExists(save)) {
        Polar old;
        bool use = ReadPolarFile(save, old);
        if (use && (old.re != polar.re || old.mach != polar.mach ||
                    old.ncrit_top != polar.ncrit_top)) {
            auto row = [this](const char* label, double a, double b,
                              bool integer) {
                out_ << ' ' << A(label, 6) << ":   ";
                if (integer)
                    out_ << I(static_cast<long>(a), 7) << std::string(25, ' ')
                         << I(static_cast<long>(b), 7);
                else
                    out_ << F(a, 12, 4) << std::string(20, ' ') << F(b, 12, 4);
                out_ << '\n';
            };
            char re_now[32];
            char re_old[32];
            std::snprintf(re_now, sizeof(re_now), "%#12.0f", polar.re);
            std::snprintf(re_old, sizeof(re_old), "%#12.0f", old.re);
            out_ << "\n               Current                         Save "
                    "file\n"
                 << "           ------------------              "
                    "------------------\n"
                 << " name  :   " << A(polar.name, 48) << A(old.name, 48)
                 << '\n'
                 << " Re    :   " << re_now << std::string(20, ' ') << re_old
                 << '\n';
            row("Mach ", polar.mach, old.mach, false);
            row("Retyp", 1, old.re_type, true);
            row("Matyp", 1, old.mach_type, true);
            row("NcritT", polar.ncrit_top, old.ncrit_top, false);
            row("NcritB", polar.ncrit_bot, old.ncrit_bot, false);
            row("xtr T", 1.0, old.xtrip_top, false);
            row("xtr B", 1.0, old.xtrip_bot, false);
            out_ << "\n Current parameters different from old save file "
                    "values.\n";
            use = AskYesNo("Set current parameters to old save file values ?");
            if (use) {
                re_ = old.re;
                viscous_ = re_ > 0.0;
                mach_ = old.mach;
                bl_.ncrit = old.ncrit_top;
                polar.re = old.re;
                polar.mach = old.mach;
                polar.ncrit_top = polar.ncrit_bot = old.ncrit_top;
            }
        }
        if (use) {
            out_ << "\n Old polar save file read in ...\n"
                 << FormatPolar(old, kProgram, version_)
                 << "\n Old polar save file available for appending\n";
            polar.points = old.points;
            polar.save_file = save;
        } else {
            out_ << "\n Polar save file will NOT be written\n";
        }
    } else {
        polar.save_file = save;
        WritePolarFile(polar);
        out_ << "\n New polar save file available\n";
    }

    const std::string dump =
        AskString("Enter  polar dump filename  OR  <return> for no file");
    if (!dump.empty())
        out_ << "\n Polar dump files are not available in Scirocco\n";
    out_ << "\n Polar dump file will NOT be written\n";

    polars_.push_back(polar);
    active_ = static_cast<int>(polars_.size()) - 1;
    pacc_ = true;
    out_ << "\n Polar accumulation enabled\n";
}

void Session::AddPolarPoint() {
    if (active_ < 0) return;
    Polar& polar = polars_[active_];
    const FoilState& state = foil_->GetState();
    PolarPoint point;
    point.alpha = state.oper.alpha;
    point.cl = state.post.cl;
    point.cm = state.post.cm;
    point.cdp = state.post.cdpi;
    if (viscous_) {
        point.cd = state.post.cd;
        const Transition top = FindTransition(state, 1);
        const Transition bottom = FindTransition(state, 0);
        point.xtr_top = top.x;
        point.xtr_bot = bottom.x;
        point.itr_top = top.node;
        point.itr_bot = bottom.node;
    }
    polar.points.push_back(point);
    out_ << "\n Point added to stored polar" << I(active_ + 1, 3) << '\n';
    if (!polar.save_file.empty()) {
        WritePolarFile(polar);
        out_ << " Point written to save file  " << A(polar.save_file, 48)
             << '\n';
    } else {
        out_ << " Save file unspecified or not available\n";
    }
    out_ << " Dump file unspecified or not available\n";
}

void Session::ListPolars(int index) const {
    if (polars_.empty()) {
        out_ << " No polars are stored\n";
        return;
    }
    out_ << '\n';
    for (std::size_t k = 0; k < polars_.size(); ++k) {
        if (index > 0 && static_cast<int>(k) + 1 != index) continue;
        out_ << kPolarRule << " Polar" << I(static_cast<long>(k) + 1, 3) << '\n'
             << FormatPolar(polars_[k], kProgram, version_);
    }
}

void Session::PolarSummary() const {
    if (polars_.empty()) {
        out_ << " No polars are stored\n";
        return;
    }
    out_ << "\n        airfoil                    Re           Mach       "
            "NcritT  NcritB  XtripT  XtripB       file\n"
         << "       ------------------------  ------------  ----------  "
            "------  ------  ------  ------    -------------------\n";
    for (std::size_t k = 0; k < polars_.size(); ++k) {
        const Polar& p = polars_[k];
        char re[32];
        std::snprintf(re, sizeof(re), "%.3fe6", p.re / 1.0e6);
        char mach[32];
        std::snprintf(mach, sizeof(mach), "%.3f", p.mach);
        out_ << I(static_cast<long>(k) + 1, 5) << "  " << A(p.name, 24) << "  "
             << A(re, 12) << "  " << A(mach, 10) << "  " << F(p.ncrit_top, 6, 2)
             << "  " << F(p.ncrit_bot, 6, 2) << "  " << F(p.xtrip_top, 6, 3)
             << "  " << F(p.xtrip_bot, 6, 3) << "    " << p.save_file << '\n';
    }
}

void Session::WritePolar(const std::vector<double>& numbers) {
    if (polars_.empty()) {
        out_ << " No polars are stored\n";
        return;
    }
    int index = 1;
    if (!numbers.empty()) {
        index = static_cast<int>(numbers[0]);
    } else if (polars_.size() > 1) {
        PolarSummary();
        AskInt("Enter index of polar to write", index);
    }
    if (index < 1 || index > static_cast<int>(polars_.size())) return;
    Polar& polar = polars_[index - 1];
    if (numbers.empty() && polars_.size() == 1) PolarSummary();
    std::string path =
        AskString("Enter polar output filename [" + polar.save_file + "]");
    if (path.empty()) path = polar.save_file;
    if (path.empty()) return;
    if (!ConfirmOverwrite(path, " Polar not written.")) return;
    std::ofstream stream(path);
    stream << FormatPolar(polar, kProgram, version_);
}

void Session::DeletePolar(const std::vector<double>& numbers) {
    if (polars_.empty()) {
        out_ << " No polars are stored\n";
        return;
    }
    int index = numbers.empty() ? 0 : static_cast<int>(numbers[0]);
    if (numbers.empty()) {
        PolarSummary();
        AskInt("Enter index of polar to delete", index);
    }
    if (index < 1 || index > static_cast<int>(polars_.size())) return;
    polars_.erase(polars_.begin() + (index - 1));
    if (active_ == index - 1) {
        active_ = -1;
        pacc_ = false;
    } else if (active_ > index - 1) {
        --active_;
    }
}

void Session::ReadPolar(const std::string& arg) {
    std::string path = Trim(arg);
    if (path.empty()) path = AskString("Enter polar filename");
    if (path.empty()) return;
    Polar polar;
    if (!ReadPolarFile(path, polar)) {
        out_ << "\n File OPEN error.  Nonexistent file:  " << path << "  \n";
        return;
    }
    polars_.push_back(polar);
    ListPolars(static_cast<int>(polars_.size()));
}

/* ------------------------------------------------------------------ */
/* Output files                                                         */
/* ------------------------------------------------------------------ */

void Session::WriteCp(const std::string& arg) {
    if (!foil_ || foil_->GetCp().size() == 0) {
        out_ << " ***  No solution available  ***\n";
        return;
    }
    std::string path = Trim(arg);
    if (path.empty()) path = AskString("Enter filename");
    if (path.empty()) return;
    const MatrixXd& x = foil_->GetFoilCoords();
    const auto& cp = foil_->GetCp();
    std::ofstream stream(path);
    stream << "#      x          Cp  \n";
    for (Eigen::Index i = x.cols() - 1; i >= 0; --i)
        stream << F(x(0, i), 12, 5) << F(cp(i), 11, 5) << '\n';
}

void Session::WriteDump(const std::string& arg) {
    if (!foil_ || foil_->GetCp().size() == 0) {
        out_ << " ***  No solution available  ***\n";
        return;
    }
    std::string path = Trim(arg);
    if (path.empty()) path = AskString("Enter filename");
    if (path.empty()) return;
    const FoilState& state = foil_->GetState();
    const MatrixXd& x = state.foil.x;
    const int n = state.foil.N;
    std::ofstream stream(path);
    stream << "#    s        x        y     Ue/Vinf    Dstar     Theta      "
              "Cf       H       H*        P         m          K          "
              "tau         Di\n";

    /* Arc length from the upper trailing edge, as XFOIL counts it */
    std::vector<double> s(n, 0.0);
    for (int i = n - 2; i >= 0; --i)
        s[i] = s[i + 1] + (x.col(i + 1) - x.col(i)).norm();

    const bool viscous = state.oper.viscous && state.glob.U.cols() >= n;
    std::vector<int> side(state.glob.U.cols(), 1);
    std::vector<scirocco::Param> params(3);
    if (viscous) {
        for (int j : state.vsol.Is[0]) side[j] = -1;
        for (int si = 0; si < 3; ++si)
            params[si] = scirocco::BuildParam(state, si);
    }
    auto surface_of = [&state](int j) {
        for (int si = 0; si < 3; ++si)
            for (int k : state.vsol.Is[si])
                if (k == j) return si;
        return 1;
    };

    for (int i = n - 1; i >= 0; --i) {
        stream << F(s[i], 10, 5) << F(x(0, i), 9, 5) << F(x(1, i), 9, 5);
        if (!viscous) {
            /* Inviscid: surface speed only */
            const double ue =
                state.isol.gam.size() > i ? -state.isol.gam(i) : 0.0;
            stream << F(ue, 9, 5) << F(0.0, 10, 6) << F(0.0, 10, 6)
                   << F(0.0, 10, 6) << F(0.0, 10, 4) << F(0.0, 10, 4)
                   << F(0.0, 9, 5) << F(0.0, 9, 5) << F(0.0, 9, 5) << '\n';
            continue;
        }
        const double ue = side[i] * state.post.ue(i);
        const double ds = state.post.ds(i);
        const double th = state.post.th(i);
        const int si = surface_of(i);
        scirocco::Param param = params[si];
        scirocco::StationParam(state, param, i);
        const Vector4d Ui = state.glob.U.col(i);
        auto [hs, hs_U] = scirocco::GetHs(Ui, param);
        stream << F(ue, 9, 5) << F(ds, 10, 6) << F(th, 10, 6)
               << F(state.post.cf(i), 10, 6) << F(th > 0 ? ds / th : 0.0, 10, 4)
               << F(hs, 10, 4) << F(ue * ue * th, 9, 5) << F(ue * ds, 9, 5)
               << F(ue * ue * ue * hs * th, 9, 5) << '\n';
    }
    if (!viscous) return;

    /* Wake: arc length continues from the trailing edge */
    double sw = s[0];
    Vector2d previous = 0.5 * (x.col(0) + x.col(n - 1));
    for (int k = 0; k < state.wake.N; ++k) {
        const int j = n + k;
        const Vector2d xw = state.wake.x.col(k);
        sw += (xw - previous).norm();
        previous = xw;
        const double ds = state.post.ds(j);
        const double th = state.post.th(j);
        stream << F(sw, 10, 5) << F(xw(0), 9, 5) << F(xw(1), 9, 5)
               << F(state.post.ue(j), 9, 5) << F(ds, 10, 6) << F(th, 10, 6)
               << F(0.0, 10, 6) << F(th > 0 ? ds / th : 0.0, 10, 4) << '\n';
    }
}

void Session::MinimumCp() const {
    if (!foil_ || foil_->GetCp().size() == 0) {
        out_ << " ***  No solution available  ***\n";
        return;
    }
    const FoilState& state = foil_->GetState();
    const MatrixXd& x = state.foil.x;
    auto minimum = [&](const scirocco::VectorXd& cp, double& xmin) {
        double cmin = 1e300;
        for (int i = 0; i < state.foil.N && i < cp.size(); ++i)
            if (cp(i) < cmin) {
                cmin = cp(i);
                xmin = x(0, i);
            }
        return cmin;
    };
    double xi = 0.0;
    const auto& cpi =
        state.post.cpi.size() > 0 ? state.post.cpi : state.post.cp;
    const double cmin_i = minimum(cpi, xi);
    out_ << "  Minimum Inviscid Cp =" << F(cmin_i, 8, 4)
         << "   at x =" << F(xi, 8, 4) << '\n';
    if (state.oper.viscous) {
        double xv = 0.0;
        const double cmin_v = minimum(state.post.cp, xv);
        out_ << "  Minimum Viscous  Cp =" << F(cmin_v, 8, 4)
             << "   at x =" << F(xv, 8, 4) << '\n';
    }
}

/* ------------------------------------------------------------------ */
/* Plots                                                                */
/* ------------------------------------------------------------------ */

void Session::SetDisplay(Display* display) {
    display_ = display;
    graphics_ = display != nullptr;
}

void Session::Show(Scene scene, PlotKind kind) {
    /* Plot-object size of the menu that draws this plot, then the zoom */
    double size = oper_size_;
    if (kind == PlotKind::Bl) size = vplo_size_;
    if (kind == PlotKind::Airfoil || kind == PlotKind::Curvature)
        size = gdes_size_;
    if (size > 0.0) scene.Scale(size / 10.0);
    scene.view = zoom_;
    plot_ = std::move(scene);
    plot_kind_ = kind;
    Reshow();
}

void Session::Reshow() {
    if (display_ != nullptr && graphics_) display_->Show(plot_);
}

bool Session::HaveSolution() const {
    return foil_ && foil_->GetCp().size() > 0;
}

Solution Session::CurrentSolution() const {
    const FoilState& state = foil_->GetState();
    double xtr_top = 1.0;
    double xtr_bot = 1.0;
    if (state.oper.viscous && state.glob.U.cols() > state.foil.N) {
        xtr_top = FindTransition(state, 1).x;
        xtr_bot = FindTransition(state, 0).x;
    }
    Solution solution = Snapshot(state, xtr_top, xtr_bot);
    solution.name = name_;
    return solution;
}

void Session::PlotCp() {
    /* Reference data is asked for each time, as XFOIL does */
    if (cp_reference_on_) ReadCpReference();
    if (force_reference_on_) ReadForceReference();
    cp_options_.version = version_;
    cp_options_.cp_reference =
        cp_reference_on_ && cp_reference_ ? &*cp_reference_ : nullptr;
    cp_options_.force_reference =
        force_reference_on_ && force_reference_ ? &*force_reference_ : nullptr;
    Show(CpPlot(CurrentSolution(), cp_options_), PlotKind::Cp);
    cp_options_.cp_reference = nullptr;
    cp_options_.force_reference = nullptr;
}

void Session::PlotPoint() {
    if (!HaveSolution()) return;
    if (!in_sequence_) {
        PlotCp();
        return;
    }
    sequence_.push_back(CurrentSolution());
    if (sequence_polar_ && pacc_ && !polars_.empty()) {
        PlotPolars();
        return;
    }
    cp_options_.version = version_;
    Show(CpSequencePlot(sequence_, cp_options_), PlotKind::CpSequence);
}

void Session::PlotCpVectors() {
    Show(CpVectorPlot(CurrentSolution()), PlotKind::CpVectors);
}

void Session::PlotBl() {
    if (!HaveSolution()) return;
    const Solution solution = CurrentSolution();
    if (!solution.viscous) return;
    Show(BlPlot(solution, bl_variable_, bl_axes_), PlotKind::Bl);
}

void Session::PlotPolars() {
    Show(PolarPlot(polars_, version_, polar_axes_, reference_polars_),
         PlotKind::Polar);
}

void Session::PlotAirfoil() {
    if (buffer_.cols() < 2) return;
    std::vector<Point> contour;
    contour.reserve(static_cast<std::size_t>(buffer_.cols()));
    for (Eigen::Index i = 0; i < buffer_.cols(); ++i)
        contour.push_back({buffer_(0, i), buffer_(1, i)});
    Show(AirfoilPlot(contour, Geometry(name_, contour), airfoil_options_),
         PlotKind::Airfoil);
}

void Session::Replot() {
    switch (plot_kind_) {
        case PlotKind::None:
            break;
        case PlotKind::Cp:
            if (HaveSolution()) PlotCp();
            break;
        case PlotKind::CpSequence:
            cp_options_.version = version_;
            Show(CpSequencePlot(sequence_, cp_options_), PlotKind::CpSequence);
            break;
        case PlotKind::CpVectors:
            if (HaveSolution()) PlotCpVectors();
            break;
        case PlotKind::Bl:
            PlotBl();
            break;
        case PlotKind::Polar:
            PlotPolars();
            break;
        case PlotKind::Airfoil:
            PlotAirfoil();
            break;
        case PlotKind::PolarAirfoils:
        case PlotKind::Curvature:
            Reshow();
            break;
    }
}

void Session::Hardcopy() {
    if (plot_kind_ == PlotKind::None) {
        out_ << " ***  No current plot  ***\n";
        return;
    }
    out_ << (hardcopies_.empty()
                 ? " Writing   PostScript to file  plot.ps ...\n"
                 : " Appending PostScript to file  plot.ps ...\n");
    hardcopies_.push_back(plot_);
    std::ofstream stream("plot.ps");
    stream << PostScript(hardcopies_);
}

void Session::AskAxis(const char* name, Axis& axis,
                      const std::vector<double>& given) {
    std::vector<double> values = given;
    if (values.empty()) {
        const std::string n = name;
        out_ << " Currently, " << n << "min," << n << "max," << n
             << "del =" << F(axis.min, 11, 4) << F(axis.max, 11, 4)
             << F(axis.del, 11, 4) << "\n Enter new  " << n << "min," << n
             << "max," << n << "del: \n"
             << std::flush;
        values = Numbers(ReadLine());
    }
    Axis changed = axis;
    if (!values.empty()) changed.min = values[0];
    if (values.size() > 1) changed.max = values[1];
    if (values.size() > 2) changed.del = values[2];
    if (changed.max > changed.min && changed.del > 0.0) axis = changed;
}

void Session::Vplo() {
    PlotBl();
    struct Choice {
        const char* name;
        BlVariable variable;
    };
    static const Choice kChoices[] = {
        {"H", BlVariable::Hk},
        {"DT", BlVariable::DstarThetaTop},
        {"DB", BlVariable::DstarThetaBottom},
        {"UE", BlVariable::Ue},
        {"CF", BlVariable::Cf},
        {"CD", BlVariable::Dissipation},
        {"N", BlVariable::Amplification},
        {"CT", BlVariable::Ctau},
        {"RT", BlVariable::RetTheta},
        {"RTL", BlVariable::LogRetTheta},
    };
    for (;;) {
        const Command command = Parse(AskCommand("..VPLO"));
        if (command.name.empty()) return;
        bool chosen = false;
        for (const Choice& choice : kChoices) {
            if (command.name == choice.name) {
                bl_variable_ = choice.variable;
                bl_axes_.y_set = false;
                PlotBl();
                chosen = true;
            }
        }
        if (chosen) continue;
        if (command.name == "?") {
            out_ << kVploMenu;
        } else if (command.name == "X" || command.name == "Y") {
            const bool x = command.name == "X";
            BlAxes current = bl_axes_;
            if (HaveSolution())
                current =
                    DefaultBlAxes(CurrentSolution(), bl_variable_, bl_axes_);
            Axis axis = x ? current.x : current.y;
            AskAxis(x ? "X" : "Y", axis, command.numbers);
            if (x) {
                bl_axes_.x = axis;
                bl_axes_.x_set = true;
            } else {
                bl_axes_.y = axis;
                bl_axes_.y_set = true;
            }
            PlotBl();
        } else if (command.name == "XPLO") {
            bl_axes_.vs_s = !bl_axes_.vs_s;
            bl_axes_.x_set = false;
            PlotBl();
        } else if (command.name == "RESE") {
            bl_axes_.x_set = false;
            bl_axes_.y_set = false;
            PlotBl();
        } else if (command.name == "GRID") {
            bl_axes_.grid = !bl_axes_.grid;
            PlotBl();
        } else if (command.name == "HARD") {
            Hardcopy();
        } else if (Is(command, {"BLOW", "B"})) {
            BlBlowup();
        } else if (command.name == "SYMB") {
            bl_axes_.symbols = !bl_axes_.symbols;
            PlotBl();
        } else if (command.name == "LABE") {
            bl_axes_.labels = !bl_axes_.labels;
            PlotBl();
        } else if (command.name == "CLIP") {
            bl_axes_.clip = !bl_axes_.clip;
            PlotBl();
        } else if (command.name == "DUMP") {
            DumpBlVariable(command.arg);
        } else if (command.name == "OVER") {
            OverlayBlVariable(command.arg);
        } else if (command.name == "SIZE") {
            SizeCommand(vplo_size_, command.numbers);
        } else if (command.name == "ANNO") {
            Annotate();
        } else if (Is(command, {"Z", "U"})) {
            ZoomCommand(command.name);
        } else if (command.name == "FRPL") {
            NotAvailable(out_, command.name);
        } else {
            NotRecognized(out_, command.name);
        }
    }
}

/* ------------------------------------------------------------------ */
/* Cursor input, zoom, plot size and annotation                         */
/* ------------------------------------------------------------------ */

bool Session::CursorReady() {
    if (display_ == nullptr || !graphics_ || plot_kind_ == PlotKind::None) {
        out_ << " ***  Cursor input needs the plot window  ***\n";
        return false;
    }
    return true;
}

std::optional<CursorHit> Session::Cursor(bool done_button) {
    out_ << std::flush;
    if (display_ == nullptr) return std::nullopt;
    return display_->Cursor(done_button);
}

bool Session::MarkCorners(Point& a, Point& b, bool absolute) {
    out_ << "\n Mark off corners of blowup area\n"
         << " (2 identical points default to current area)\n";
    Point* corners[2] = {&a, &b};
    for (int k = 0; k < 2; ++k) {
        const std::optional<CursorHit> hit = Cursor(false);
        if (!hit) return false;
        *corners[k] = hit->page;
        const Point user = plot_.user.set ? plot_.user.ToUser(hit->page)
                                          : Point{hit->page.x / kUnitsPerInch,
                                                  hit->page.y / kUnitsPerInch};
        if (absolute)
            out_ << " Pt" << I(k + 1, 2) << "  usr x,y" << F(user.x, 13, 6)
                 << F(user.y, 13, 6) << "   abs X,Y"
                 << F(hit->page.x / kUnitsPerInch, 13, 6)
                 << F(hit->page.y / kUnitsPerInch, 13, 6) << '\n';
        else
            out_ << " x,y = " << ListDirected(user.x) << ' '
                 << ListDirected(user.y) << '\n';
    }
    return std::hypot(a.x - b.x, a.y - b.y) > 1e-6;
}

bool Session::ZoomCommand(const std::string& name) {
    if (name == "U") {
        zoom_ = Rect{};
        plot_.view = zoom_;
        Reshow();
        return true;
    }
    if (name != "Z") return false;
    if (!CursorReady()) return true;
    Point a;
    Point b;
    if (!MarkCorners(a, b, true)) return true;
    zoom_ = AspectRect(a, b);
    plot_.view = zoom_;
    Reshow();
    return true;
}

void Session::SizeCommand(double& size, const std::vector<double>& given) {
    double value = size;
    if (!given.empty()) {
        value = given[0];
    } else {
        out_ << " Current plot-object size = " << ListDirected(size) << '\n';
        AskReal("Enter new plot-object size", value);
    }
    if (value > 0.0) size = value;
}

void Session::Annotate() {
    if (!CursorReady()) return;
    constexpr double kText = 17.5;
    constexpr double kSymbol = 10.0;
    out_ << kAnnoMenu;
    for (;;) {
        out_ << " Select option or <return>:  " << std::flush;
        const std::string line = Trim(ReadLine());
        if (line.empty()) return;
        const char option = static_cast<char>(
            std::toupper(static_cast<unsigned char>(line[0])));
        if (option == 'C' || option == 'S' || option == 'M') {
            out_ << " Click on lower left point of character string...\n";
            const std::optional<CursorHit> hit = Cursor(false);
            if (!hit) return;
            out_ << " Enter character string:  " << std::flush;
            std::string text = ReadLine();
            out_ << '\n';
            if (option == 'M') text = MathCharacters(text);
            Label label;
            label.at = hit->page;
            label.height = kText;
            label.text = text;
            label.width = anno_width_;
            label.slant = option == 'S';
            plot_.labels.push_back(label);
            Reshow();
        } else if (option == 'P') {
            out_ << " Click on symbol locations ...\n";
            for (;;) {
                const std::optional<CursorHit> hit = Cursor(true);
                if (!hit || hit->key == 'D') break;
                plot_.Symbol(hit->page.x, hit->page.y, kSymbol, anno_symbol_);
                Reshow();
            }
            out_ << '\n';
        } else if (option == 'L' || option == 'A') {
            out_ << " Click on line points, twice on last point...\n";
            std::vector<Point> points;
            for (;;) {
                const std::optional<CursorHit> hit = Cursor(true);
                if (!hit || hit->key == 'D') break;
                if (!points.empty() &&
                    std::hypot(hit->page.x - points.back().x,
                               hit->page.y - points.back().y) < 3.0)
                    break;
                points.push_back(hit->page);
            }
            if (points.size() >= 2) {
                plot_.Line(points);
                if (option == 'A') {
                    const Point& tip = points.back();
                    const Point& from = points[points.size() - 2];
                    const double length =
                        std::hypot(tip.x - from.x, tip.y - from.y);
                    const double dx = (tip.x - from.x) / length;
                    const double dy = (tip.y - from.y) / length;
                    constexpr double kHead = 16.0;
                    constexpr double kHalf = 4.5;
                    const Point base{tip.x - dx * kHead, tip.y - dy * kHead};
                    plot_.Line({{base.x - dy * kHalf, base.y + dx * kHalf},
                                tip,
                                {base.x + dy * kHalf, base.y - dx * kHalf}});
                }
                Reshow();
            }
            out_ << '\n';
        } else if (option == 'W') {
            out_ << " Enter new character width factor (currently ="
                 << F(anno_width_, 7, 3) << "):  " << std::flush;
            const std::vector<double> v = Numbers(ReadLine());
            out_ << '\n';
            if (!v.empty() && v[0] > 0.0) anno_width_ = v[0];
        } else if (option == 'T') {
            out_ << '\n'
                 << kSymbolTable << "\n Enter new symbol type (currently ="
                 << I(anno_symbol_, 2) << "):  " << std::flush;
            const std::vector<double> v = Numbers(ReadLine());
            out_ << '\n';
            if (!v.empty() && v[0] >= 0.0 && v[0] <= 13.0)
                anno_symbol_ = static_cast<int>(v[0]);
        } else {
            out_ << kAnnoMenu;
        }
    }
}

/* ------------------------------------------------------------------ */
/* Boundary-layer velocity profiles                                     */
/* ------------------------------------------------------------------ */

namespace {

void PrintProfile(std::ostream& out, const BlProfile& profile) {
    out << " x y =" << F(profile.wall.x, 8, 4) << F(profile.wall.y, 8, 4)
        << "    Delta* =" << E(profile.ds, 12, 4)
        << "    Rtheta =" << F(profile.ret, 10, 2)
        << "    Hk =" << F(profile.hk, 9, 4) << '\n';
}

}  // namespace

void Session::Profiles(const std::vector<double>& given) {
    if (!HaveSolution()) return;
    const Solution solution = CurrentSolution();
    if (!solution.viscous) return;
    int count = 21;
    if (!given.empty() && given[0] >= 2.0)
        count = static_cast<int>(given[0]);
    else
        out_ << " Using default number of profiles:" << I(count, 12) << '\n';
    if (plot_kind_ != PlotKind::Cp) PlotCp();
    const double arc = ArcLength(solution);
    for (int k = 0; k < count; ++k) {
        BlProfile profile;
        if (!ProfileAtArc(solution, arc * k / (count - 1), profile_weight_,
                          profile))
            continue;
        PrintProfile(out_, profile);
        DrawProfile(plot_, profile);
    }
    Reshow();
}

void Session::ProfilesAtCursor() {
    if (!HaveSolution()) return;
    const Solution solution = CurrentSolution();
    if (!solution.viscous) return;
    if (plot_kind_ != PlotKind::Cp) PlotCp();
    if (!CursorReady()) return;
    out_ << "  \n Locate profiles with cursor, type \"D\" when done...\n";
    for (;;) {
        const std::optional<CursorHit> hit = Cursor(true);
        if (!hit || hit->key == 'D') break;
        BlProfile profile;
        if (!ProfileNear(solution, plot_.foil.ToUser(hit->page),
                         profile_weight_, profile))
            continue;
        PrintProfile(out_, profile);
        DrawProfile(plot_, profile);
        Reshow();
    }
}

void Session::ProfileWeight(const std::vector<double>& given) {
    double weight = profile_weight_;
    if (!given.empty()) {
        weight = given[0];
    } else {
        out_ << "\n Current u/Qinf profile plot weight =" << E(weight, 12, 4)
             << '\n';
        AskReal("Enter new plot weight", weight);
    }
    if (weight > 0.0) profile_weight_ = weight;
}

/* ------------------------------------------------------------------ */
/* Reference data                                                       */
/* ------------------------------------------------------------------ */

void Session::ReadCpReference() {
    out_ << "\n Enter Cp vs x data filename " << cp_reference_file_ << '\n'
         << std::flush;
    const std::string answer = Trim(ReadLine());
    if (!answer.empty()) cp_reference_file_ = answer;
    std::vector<std::vector<Point>> blocks;
    if (cp_reference_file_.empty() ||
        !ReadPairBlocks(cp_reference_file_, blocks)) {
        out_ << " GETXYL: File OPEN error.\n";
        cp_reference_.reset();
        return;
    }
    CpReference reference;
    for (const auto& block : blocks)
        reference.points.insert(reference.points.end(), block.begin(),
                                block.end());
    reference.label = AskString("Enter data description label");
    cp_reference_ = reference;
}

void Session::ReadForceReference() {
    out_ << " Enter reference  Mach, Re, Alpha, CL, CD, CM:\n" << std::flush;
    const std::vector<double> v = Numbers(ReadLine());
    ForceReference reference =
        force_reference_ ? *force_reference_ : ForceReference{};
    double* fields[6] = {&reference.mach, &reference.re, &reference.alpha,
                         &reference.cl,   &reference.cd, &reference.cm};
    for (std::size_t k = 0; k < v.size() && k < 6; ++k) *fields[k] = v[k];
    reference.label = AskString("Enter reference force data label");
    force_reference_ = reference;
}

void Session::ReadReferencePolar(const std::string& arg) {
    std::string path = Trim(arg);
    if (path.empty()) path = AskString("Enter reference polar filename");
    if (path.empty()) return;
    std::vector<std::vector<Point>> blocks;
    if (!ReadPairBlocks(path, blocks)) {
        out_ << "\n File OPEN error.  Nonexistent file:  " << path << "  \n";
        return;
    }
    ReferencePolar reference;
    std::vector<Point>* sets[4] = {&reference.cd_cl, &reference.alpha_cl,
                                   &reference.alpha_cm, &reference.xtr_cl};
    for (std::size_t k = 0; k < blocks.size() && k < 4; ++k)
        *sets[k] = blocks[k];
    reference.label = AskString("Enter label for reference polar");
    reference_polars_.push_back(reference);
    if (plot_kind_ == PlotKind::Polar) PlotPolars();
}

void Session::DeleteReferencePolar(const std::vector<double>& given) {
    if (reference_polars_.empty()) {
        out_ << " No reference polars are stored\n";
        return;
    }
    int which = -1;
    if (!given.empty()) {
        which = static_cast<int>(given[0]);
    } else {
        out_ << "\n        reference polar                          \n"
             << "       ------------------------------------------\n";
        for (std::size_t k = 0; k < reference_polars_.size(); ++k)
            out_ << I(static_cast<long>(k) + 1, 5) << "  "
                 << A(reference_polars_[k].label, 48) << '\n';
        AskInt("Specify ref. polar to delete (0 = all, -1 = abort)", which);
    }
    if (which == 0)
        reference_polars_.clear();
    else if (which > 0 &&
             static_cast<std::size_t>(which) <= reference_polars_.size())
        reference_polars_.erase(reference_polars_.begin() + which - 1);
    if (plot_kind_ == PlotKind::Polar) PlotPolars();
}

void Session::PolarAxesCommand() {
    PolarAxes axes = DefaultPolarAxes(polars_, polar_axes_);
    out_
        << " Enter new axis annotations, or <return> to leave unchanged...\n\n";
    struct Row {
        const char* label;
        Axis* axis;
    };
    const Row rows[] = {{"   Alpha  min, max, delta:", &axes.alpha},
                        {"     CL   min, max, delta:", &axes.cl},
                        {"     CD   min, max, delta:", &axes.cd},
                        {"    -CM   min, max, delta:", &axes.cm}};
    for (const Row& row : rows) {
        out_ << row.label << F(row.axis->min, 11, 5) << F(row.axis->max, 11, 5)
             << F(row.axis->del, 11, 5) << '\n'
             << std::flush;
        const std::vector<double> v = Numbers(ReadLine());
        Axis changed = *row.axis;
        if (!v.empty()) changed.min = v[0];
        if (v.size() > 1) changed.max = v[1];
        if (v.size() > 2) changed.del = v[2];
        if (changed.max > changed.min && changed.del > 0.0) *row.axis = changed;
    }
    axes.set = true;
    polar_axes_ = axes;
    if (plot_kind_ == PlotKind::Polar) PlotPolars();
}

/* ------------------------------------------------------------------ */
/* VPLO extras                                                          */
/* ------------------------------------------------------------------ */

void Session::BlBlowup() {
    if (plot_kind_ != PlotKind::Bl || !CursorReady()) return;
    Point a;
    Point b;
    if (!MarkCorners(a, b, false)) return;
    const Point ua = plot_.user.ToUser(a);
    const Point ub = plot_.user.ToUser(b);
    bl_axes_.x = NiceAxisFor(std::min(ua.x, ub.x), std::max(ua.x, ub.x), 7);
    bl_axes_.y = NiceAxisFor(std::min(ua.y, ub.y), std::max(ua.y, ub.y), 8);
    bl_axes_.x_set = true;
    bl_axes_.y_set = true;
    PlotBl();
}

void Session::DumpBlVariable(const std::string& arg) {
    if (!HaveSolution()) return;
    const Solution solution = CurrentSolution();
    if (!solution.viscous) return;
    std::string path = Trim(arg);
    if (path.empty()) path = AskString("Enter filename");
    if (path.empty()) return;
    std::ofstream stream(path);
    stream << "#" << (bl_axes_.vs_s ? "     s" : "     x")
           << "          value\n";
    for (const auto& curve : BlData(solution, bl_variable_, bl_axes_.vs_s)) {
        for (const Point& p : curve)
            stream << F(p.x, 10, 6) << ' ' << E(p.y, 14, 6) << '\n';
        stream << '\n';
    }
}

void Session::OverlayBlVariable(const std::string& arg) {
    std::string path = Trim(arg);
    if (path.empty()) path = AskString("Enter filename");
    if (path.empty()) return;
    std::vector<std::vector<Point>> blocks;
    if (!ReadPairBlocks(path, blocks)) {
        out_ << "\n File OPEN error.  Nonexistent file:  " << path << "  \n";
        return;
    }
    bl_axes_.overlay.clear();
    for (const auto& block : blocks)
        bl_axes_.overlay.insert(bl_axes_.overlay.end(), block.begin(),
                                block.end());
    PlotBl();
}

/* ------------------------------------------------------------------ */
/* GDES plot commands                                                   */
/* ------------------------------------------------------------------ */

namespace {

std::vector<Point> Contour(const MatrixXd& x) {
    std::vector<Point> contour;
    contour.reserve(static_cast<std::size_t>(x.cols()));
    for (Eigen::Index i = 0; i < x.cols(); ++i)
        contour.push_back({x(0, i), x(1, i)});
    return contour;
}

}  // namespace

void Session::GdesBlowup() {
    if (plot_kind_ != PlotKind::Airfoil || !CursorReady()) return;
    Point a;
    Point b;
    if (!MarkCorners(a, b, false)) return;
    const Point ua = plot_.user.ToUser(a);
    const Point ub = plot_.user.ToUser(b);
    airfoil_options_.window = {std::min(ua.x, ub.x), std::min(ua.y, ub.y),
                               std::max(ua.x, ub.x), std::max(ua.y, ub.y)};
    airfoil_options_.window_set = true;
    airfoil_options_.inches = false;
    PlotAirfoil();
}

void Session::GdesWindow() {
    if (plot_kind_ != PlotKind::Airfoil || !CursorReady()) return;
    out_ << "  \n Type I,O,P to In,Out,Pan with cursor...\n";
    Rect w = airfoil_options_.window_set ? airfoil_options_.window
                                         : AirfoilWindow(Contour(buffer_));
    for (;;) {
        const std::optional<CursorHit> hit = Cursor(true);
        if (!hit) break;
        const char key = hit->key;
        if (key != 0 && key != 'I' && key != 'O' && key != 'P') break;
        const Point u = plot_.user.ToUser(hit->page);
        const double factor = key == 'I' ? 0.5 : key == 'O' ? 2.0 : 1.0;
        const double hw = 0.5 * (w.x1 - w.x0) * factor;
        const double hh = 0.5 * (w.y1 - w.y0) * factor;
        w = {u.x - hw, u.y - hh, u.x + hw, u.y + hh};
        airfoil_options_.window = w;
        airfoil_options_.window_set = true;
        airfoil_options_.inches = false;
        PlotAirfoil();
    }
}

void Session::GdesOverlay(const std::string& arg) {
    std::string path = Trim(arg);
    if (path.empty()) path = AskString("Enter filename");
    if (path.empty()) return;
    AirfoilFile file;
    std::string error;
    if (!ReadAirfoilFile(path, file, error)) {
        out_ << "\n File OPEN error.  Nonexistent file:  " << path << "  \n";
        return;
    }
    std::string name = file.name;
    if (file.labeled) {
        out_ << "\n Labeled airfoil file.  Name:  " << A(name, 32) << '\n';
    } else {
        out_ << "\n Plain airfoil file\n";
        name = AskString("Enter airfoil name");
    }
    const MatrixXd points =
        IsCounterclockwise(file.x) ? Reversed(file.x) : file.x;
    PrintThicknessCamber(points);
    airfoil_options_.overlay = Contour(points);
    airfoil_options_.overlay_info = Geometry(name, airfoil_options_.overlay);
    PlotAirfoil();
}

void Session::Curvatures() {
    if (buffer_.cols() < 3) return;
    const std::vector<Point> contour = Contour(buffer_);
    const std::vector<double> k = Curvature(contour);
    const std::size_t n = contour.size();
    /* Arc length from the upper trailing edge, XFOIL's first point */
    std::vector<double> s(n, 0.0);
    for (std::size_t i = n - 1; i-- > 0;)
        s[i] = s[i + 1] + std::hypot(contour[i].x - contour[i + 1].x,
                                     contour[i].y - contour[i + 1].y);
    std::size_t le = 0;
    std::size_t kmax = 0;
    std::size_t kmin = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (contour[i].x < contour[le].x) le = i;
        if (k[i] > k[kmax]) kmax = i;
        if (k[i] < k[kmin]) kmin = i;
    }
    auto stat = [&](const char* name, std::size_t i) {
        out_ << name << ListDirected(k[i]) << "  at S =  " << ListDirected(s[i])
             << "  at X =  " << ListDirected(contour[i].x) << '\n';
    };
    double exponent = 1.0;
    for (;;) {
        out_ << "  \n Airfoil curvature (yellow-upper, cyan-lower) \n";
        stat("  CVLE  =  ", le);
        stat("  CVmax =  ", kmax);
        stat("  CVmin =  ", kmin);
        out_ << "  \n Enter C for curvature    plot\n"
             << " Enter N for curvature**N plot\n"
             << " Hit <cr> to exit\n";
        const Command command = Parse(AskCommand("..CPLO"));
        if (command.name.empty()) {
            if (plot_kind_ == PlotKind::Curvature) PlotAirfoil();
            return;
        }
        if (command.name == "C") {
            Show(CurvaturePlot(contour, 1.0), PlotKind::Curvature);
        } else if (command.name == "N") {
            if (!command.numbers.empty())
                exponent = command.numbers[0];
            else
                AskReal("Enter curvature exponent N", exponent);
            Show(CurvaturePlot(contour, exponent), PlotKind::Curvature);
        } else if (Is(command, {"Z", "U"})) {
            ZoomCommand(command.name);
        }
    }
}

void Session::Distance() {
    if (!CursorReady()) return;
    out_ << "\n Click mouse or hit a key on each point\n\n";
    Point u[2];
    for (Point& p : u) {
        const std::optional<CursorHit> hit = Cursor(false);
        if (!hit) return;
        p = plot_.user.set ? plot_.user.ToUser(hit->page)
                           : Point{hit->page.x / kUnitsPerInch,
                                   hit->page.y / kUnitsPerInch};
    }
    const double dx = u[1].x - u[0].x;
    const double dy = u[1].y - u[0].y;
    out_ << " x1 =" << F(u[0].x, 10, 6) << "    y1 =" << F(u[0].y, 10, 6)
         << "\n x2 =" << F(u[1].x, 10, 6) << "    y2 =" << F(u[1].y, 10, 6)
         << "\n dx =" << F(dx, 10, 6) << "    dy =" << F(dy, 10, 6)
         << "    ds =" << F(std::hypot(dx, dy), 10, 6) << '\n';
}

void Session::TickSize(const std::vector<double>& given) {
    double size = airfoil_options_.tick_size;
    if (!given.empty()) {
        size = given[0];
    } else {
        out_ << " Current tick-mark size (as fraction of perimeter) = "
             << ListDirected(size) << '\n';
        AskReal("Enter new tick-mark size", size);
    }
    if (size > 0.0) airfoil_options_.tick_size = size;
    if (airfoil_options_.ticks) PlotAirfoil();
}

}  // namespace fauxfoil
