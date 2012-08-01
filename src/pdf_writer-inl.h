/*
  Copyright 2010-2012 Andreas Biegert, Christof Angermueller

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CS_PDF_WRITER_INL_H_
#define CS_PDF_WRITER_INL_H_

#include "pdf_writer.h"

namespace cs {

template<class Abc>
AlignmentPdfWriter<Abc>::AlignmentPdfWriter(const Alignment<Abc>& ali)
        : ali_(ali) {
    // Generate temporary filename for FASTA alignment
    char tempfilename[] = "/tmp/alignment_XXXXXX";
    const int captured_fd = mkstemp(tempfilename);
    if (!captured_fd) throw Exception("Unable to create unique alignment filename!");
    tempfile_ = tempfilename;
}

template<class Abc>
void AlignmentPdfWriter<Abc>::WriteBody(FILE* fp) {
    FILE* fali = fopen(tempfile_.c_str(), "w");
    ali_.Write(fali, FASTA_ALIGNMENT);
    fclose(fali);

    fprintf(fp, "\\begin{texshade}{%s}\n", tempfile_.c_str());
    fputs("\\shadingmode{functional}\n", fp);
    fputs("\\funcgroup{Adenine}{A}{Black}{dnacolA!80!white}{upper}{bf}\n", fp);
    fputs("\\funcgroup{Cytosine}{C}{Black}{dnacolC!80!white}{upper}{bf}\n", fp);
    fputs("\\funcgroup{Guanine}{G}{Black}{dnacolG!80!white}{upper}{bf}\n", fp);
    fputs("\\funcgroup{Thymine}{T}{Black}{dnacolT!80!white}{upper}{bf}\n", fp);
    fputs("\\gapcolors{black}{nearlywhite}\n", fp);
    fputs("\\setsize{residues}{normalsize}\n", fp);
    fputs("\\setfamily{residues}{tt}\n", fp);
    fputs("\\setseries{residues}{bf}\n", fp);
    fputs("\\setshape{residues}{up}\n", fp);
    fprintf(fp, "\\residuesperline*{%zu}\n", ali_.ncols());
    fputs("\\seqtype{N}\n", fp);
    fputs("\\shadeallresidues\n", fp);
    fputs("\\gapchar{-}\n", fp);
    fputs("\\setfamily{features}{tt}\n", fp);
    fputs("\\setshape{features}{up}\n", fp);
    fputs("\\hidenumbering\n", fp);
    fputs("\\hidenames\n", fp);
    fputs("\\end{texshade}\n", fp);
}


template<class Abc>
ProfilePdfWriter<Abc>::ProfilePdfWriter(const Profile<Abc>& profile)
        : col_height(10.0),
          hide_numbering(false),
          show_chars(true),
          pmin(0.005),
          profile_(profile) {
    Normalize(profile_, 1.0);
}

template<class Abc>
void ProfilePdfWriter<Abc>::WriteBody(FILE* fp) {

    fputs("\\sffamily\n", fp);
    fputs("\\begin{tikzpicture}[x=2.5ex,y=2.5ex]\n", fp);
    fputs("\\tikzstyle{sRect}=[rectangle,draw=black,semithick];\n", fp);
    fputs("\\tikzstyle{sChar}=[inner sep=0];\n", fp);
    fputs("\\tikzstyle{sNumbering}=[rotate=90,anchor=east];\n", fp);
    fputs("\\tikzstyle{sLogo}=[rectangle,draw=black,inner sep=0pt,minimum size=0mm,"
          "anchor=south,semithick];\n", fp);
    fputs("\\definecolor{darkgray}{rgb}{0.3,0.3,0.3};\n", fp);


    for (size_t i = 0; i < profile_.length(); ++i) {
        // Sort residues by functional group and probability
        std::vector<GroupProbResTuple> tuples;
        for (size_t a = 0; a < Abc::kSize; ++a)
            if (profile_[i][a] > pmin)
                tuples.push_back(boost::make_tuple(-Abc::kFuncGroup[a], profile_[i][a], a));
        std::sort(tuples.begin(), tuples.end());
        std::reverse(tuples.begin(), tuples.end());

        // Draw profile logos
        double x = i;
        double y = 0;
        if (!hide_numbering) {
            if (profile_.length() < 20 && profile_.length() % 2 == 1) {
                int wlen = (profile_.length() - 1) / 2;
                int idx = abs(static_cast<int>(i) - wlen);
                fprintf(fp, "\\node[anchor=north] at (%.4f,%.4f) {\\sffamily \\textcolor{%s}{%d}};\n",
                        x + 0.5, y, idx == 0 ? "red" : "black", idx);
            } else {
                fprintf(fp, "\\node[sNumbering] at (%.4f,%.4f) {\\sffamily %zu};\n", x + 0.5, y, i + 1);
            }
        }
        if (show_chars) {
            size_t max = 0;
            for (size_t a = 1; a < Abc::kSize; ++a)
                if (profile_[i][a] > profile_[i][max])
                    max = a;
            fprintf(fp, "\\filldraw [sRect,fill=darkgray] (%.4f,%.4f) rectangle +(1,1);\n", x, y);
            fprintf(fp, "\\node[sChar] at (%.4f, %.4f) {\\bf \\sffamily \\textcolor{white}{%c}};\n", 
                    x + 0.5, y + 0.5, Abc::kIntToChar[max]);
            y += 1.0;
        }
        for (size_t k = 0; k < tuples.size(); ++k) {
            double prob  = boost::get<1>(tuples[k]);
            uint8_t res  = boost::get<2>(tuples[k]);
            double height = prob * col_height;
            fprintf(fp, "\\filldraw[sRect,fill=%scol%c] (%.4f,%.4f) rectangle +(1,%.4f);\n",
                    Abc::kName, Abc::kIntToChar[res], x, y, height);
            if (prob > 0.08)
                fprintf(fp, "\\node at (%.4f,%.4f) {\\bf \\sffamily %c};\n",
                        x + 0.5, y + 0.5 * height, Abc::kIntToChar[res]);
            y += height;
        }
    }

    fputs("\\end{tikzpicture}\n", fp);
}


template<class AS, class Abc>
StateProfilePdfWriter<AS, Abc>::StateProfilePdfWriter(const Profile<AS>& profile,
                                                      const ContextLibrary<Abc>& lib)
        : hide_numbering(false),
          pmin(0.001),
          profile_(profile),
          lib_(lib) {
    Normalize(profile_, 1.0);
}

template<class AS, class Abc>
void StateProfilePdfWriter<AS, Abc>::WriteBody(FILE* fp) {
    DefineColors(fp);

    fputs("\\sffamily\n", fp);
    fprintf(fp, "\\begin{tikzpicture}[x=%.5fcm,y=%.5fcm,\n", kResWidth, kScaleY);
    fputs("numbering/.style={rotate=90,anchor=east},\n", fp);
    fputs("logo/.style={rectangle,draw=black,inner sep=0pt,minimum size=0mm,"
          "anchor=south,semithick}]\n", fp);

    for (size_t i = 0; i < profile_.length(); ++i) {
        // Sort abstract states by color and probability
        std::vector<ColorProbStateTuple> tuples;
        for (size_t k = 0; k < AS::kSize; ++k)
            if (profile_[i][k] > pmin)
                tuples.push_back(boost::make_tuple(lib_[k].color, profile_[i][k], k));
        std::sort(tuples.begin(), tuples.end());
        std::reverse(tuples.begin(), tuples.end());

        // Draw profile logos
        double x = i - 0.5;
        double y = 0;
        if (!hide_numbering)
            fprintf(fp, "\\node[numbering] at (%zu,0) {\\tiny %zu};\n", i, i + 1);
        for (size_t k = 0; k < tuples.size(); ++k) {
            double red   = 100 * boost::get<0>(tuples[k]).red;
            double green = 100 * boost::get<0>(tuples[k]).green;
            double blue  = 100 * boost::get<0>(tuples[k]).blue;
            double prob  = boost::get<1>(tuples[k]);
            uint8_t as   = boost::get<2>(tuples[k]);
            fprintf(fp, "\\filldraw [logo,fill=col%.0f%.0f%.0f] (%.4f,%.4f) "
                    "rectangle +(%.4f,%.4f);\n", red, green, blue, x, y, 1.0, prob * 20);
            if (prob > 0.08)
                fprintf(fp, "\\node at (%zu,%.4f) {\\tiny \\bf \\sffamily %c};\n",
                        i, y + prob * 10, AS::kIntToChar[as]);
            y += prob * 20;
        }
    }

    fputs("\\end{tikzpicture}\n", fp);
}

template<class AS, class Abc>
void StateProfilePdfWriter<AS, Abc>::DefineColors(FILE* fp) const {
    std::set<Color> colors;
    for (size_t k = 0; k < lib_.size(); ++k)
        colors.insert(lib_[k].color);
    for (std::set<Color>::iterator c = colors.begin(); c != colors.end(); ++c)
        fprintf(fp, "\\definecolor{col%.0f%.0f%.0f}{rgb}{%.2f,%.2f,%.2f}\n",
                100 * c->red, 100 * c->green, 100 * c->blue,
                c->red, c->green, c->blue);
}


template<class Abc>
void ContextLibraryPdfWriter<Abc>::WriteBody(FILE* fp) {
    DefineColors(fp);

    fputs("\\sffamily\n", fp);
    fprintf(fp, "\\begin{tikzpicture}[x=%.5fcm,y=%.5fcm,\n", kResWidth, kScaleY);
    fputs("colidx/.style={anchor=north},\n", fp);
    fprintf(fp, "pname/.style={rectangle,thin,draw=black,inner sep=0.2mm,minimum "
            "size=%.5fcm,anchor=north,rounded corners=1pt},\n", kResWidth);
    fputs("logo/.style={rectangle,draw=black,inner sep=0pt,minimum size=0mm,"
          "anchor=south,semithick}]\n", fp);

    const size_t kmax = MIN(2000, static_cast<int>(lib_.size()));
    const double kNameHeight = 3.0;
    const double kRulerHeight = 4.0;
    const double kHistoHeight = 20.0;
    const double kProfileHeight = kNameHeight + kRulerHeight + kHistoHeight;
    double x = 0.0;
    double y = 0.0;

    for (size_t k = 0; k < kmax; ++k) {
        const ContextProfile<Abc>& cp = lib_[k];
        fprintf(fp, "\\message{%zu %s}\n", k, cp.name.c_str());
        if (k % ncols == 0) { x = 0.0; y -= kProfileHeight + margin; }

        // Draw name badge
        fprintf(fp, "\\node[pname,fill=col%.0f%.0f%.0f] at (%.4f,%.4f) "
                "{\\tiny \\bf \\sffamily %s};\n", 100 * cp.color.red,
                100 * cp.color.green, 100 * cp.color.blue, x + lib_.center(),
                y, cp.name.empty() ? strprintf("%zu", k).c_str() : cp.name.c_str());

        for (size_t i = 0; i < lib_.wlen(); ++i) {
            double yi = y - kNameHeight - kHistoHeight;

            // Sort residues by functional group and probability
            std::vector<GroupProbResTuple> tuples;
            for (size_t a = 0; a < Abc::kSize; ++a)
                if (cp.probs[i][a] > pmin)
                    tuples.push_back(boost::make_tuple(-Abc::kFuncGroup[a],
                                                       cp.probs[i][a], a));
            std::sort(tuples.begin(), tuples.end());
            std::reverse(tuples.begin(), tuples.end());

            // Draw index postion
            if (lib_.wlen() > 1) {
                int idx = abs(static_cast<int>(i) - static_cast<int>(lib_.center()));
                fprintf(fp, "\\node[colidx] at (%.4f,%.4f) {\\tiny "
                        "\\textcolor{%s}{%s}};\n", x, yi, idx == 0 ? "red" : "black",
                        strprintf("%d", idx).c_str());
            }

            // Draw profile logos
            for (size_t t = 0; t < tuples.size(); ++t) {
                double prob  = boost::get<1>(tuples[t]);
                uint8_t res  = boost::get<2>(tuples[t]);
                fprintf(fp, "\\filldraw[logo,fill=%scol%c] (%.4f,%.4f) rectangle "
                        "+(1,%.4f);\n",  Abc::kName, Abc::kIntToChar[res],
                        x - 0.5, yi, prob * kHistoHeight);
                if (prob > 0.08)
                    fprintf(fp, "\\node at (%.4f,%.4f) {\\tiny \\bf \\sffamily %c};\n",
                            x, yi + prob * 0.5 * kHistoHeight, Abc::kIntToChar[res]);
                yi += prob * kHistoHeight;
            }

            x += 1.0;  // advance to next profile column
        }
        x += margin;  // add horizontal margin between profiles
    }

    fputs("\\end{tikzpicture}\n", fp);
}

template<class Abc>
void ContextLibraryPdfWriter<Abc>::DefineColors(FILE* fp) const {
    std::set<Color> colors;
    for (size_t k = 0; k < lib_.size(); ++k)
        colors.insert(lib_[k].color);
    for (std::set<Color>::iterator c = colors.begin(); c != colors.end(); ++c)
        fprintf(fp, "\\definecolor{col%.0f%.0f%.0f}{rgb}{%.2f,%.2f,%.2f}\n",
                100 * c->red, 100 * c->green, 100 * c->blue,
                c->red, c->green, c->blue);
}

template<class Abc>
std::vector<size_t> CrfStatePdfWriter<Abc>::SortRes(const double* col) const {
    std::vector<GroupProbResTuple> tuples;
    for (size_t a = 0; a < Abc::kSize; ++a) {
        if (col[a] > prob_min) {
            tuples.push_back(boost::make_tuple(-Abc::kFuncGroup[a], col[a], a));
        }
    }
    std::sort(tuples.begin(), tuples.end());
    std::reverse(tuples.begin(), tuples.end());
    std::vector<size_t> res;
    for (size_t k = 0; k < tuples.size(); ++k) {
        res.push_back(boost::get<2>(tuples[k]));
    }
    return res;
}

template<class Abc>
void ProbCrfStatePdfWriter<Abc>::WriteBody(FILE* fp) {
    const CrfState<Abc>* crf_state = ProbCrfStatePdfWriter<Abc>::crf_state;
    const size_t length = crf_state->length();
    const size_t center = 0.5 * (length - 1);

    double col_weights[length];
    if (weight_center > 0.0 && weight_decay > 0.0) {
        col_weights[center] = weight_center;
        for (size_t i = 1; i <= center; ++i) {
            col_weights[center - i] = col_weights[center - i + 1] * weight_decay;
            col_weights[center + i] = col_weights[center - i];
        }
    } else {
        for (size_t i = 0; i < length; ++i) {
            const double* col = crf_state->context_weights[i];
            std::vector<double> vcol(col, col + Abc::kSize);
            sort(vcol.begin(), vcol.end());
            size_t m = 0.5 * vcol.size();
            double median = 0.5 * (vcol[m] + vcol[m - 1]);
            double mad = 0.0;
            for (size_t a = 0; a < Abc::kSize; ++a) {
                mad += fabs(col[a] - median);
            }
            mad /= Abc::kSize;
            col_weights[i] = mad;
        }
    }

    Profile<Abc> context_probs(length);
    for (size_t i = 0; i < length; ++i) {
        const double* col = crf_state->context_weights[i];
        double norm = 0.0;
        double tmp[Abc::kSize];
        for (size_t a = 0; a < Abc::kSize; ++a) {
            tmp[a] = pow(exp(col[a]), col_weights[i]);
            norm += tmp[a];
        }
        for (size_t a = 0; a < Abc::kSize; ++a) {
            context_probs[i][a] = tmp[a] / norm;
        }
    }

    ProfileColumn<Abc> pc_probs;
    double pc_sum = 0.0;
    for (size_t a = 0; a < Abc::kSize; ++a) {
        pc_sum += exp(crf_state->pc_weights[a]);
    }
    for (size_t a = 0; a < Abc::kSize; ++a) {
        pc_probs[a] = exp(crf_state->pc_weights[a]) / pc_sum;
    }

    fputs("\\begin{tikzpicture}[x=2.5ex, y=2.5ex]\n", fp);
    fputs("\\tikzstyle{sRect}=[rectangle,draw=black,semithick];\n", fp);
    fputs("\\tikzstyle{sChar}=[inner sep=0];\n", fp);
    fputs("\\definecolor{darkgray}{rgb}{0.3,0.3,0.3};\n", fp);

    // Name
    if (show_name) {
        fprintf(fp, "\\node [anchor=west] at (0,%.4f) {\\bf \\sffamily %s};\n",
                col_height + 1, crf_state->name.c_str());
    }

    double x;
    double y;

    // Context weights
    for (size_t i = 0; i < length; ++i) {
        y = 0.0;
        x = i;
        fprintf(fp, "\\filldraw [sRect, fill=%s] (%.4f,%.4f) rectangle +(1,-1);\n",
                i == center ? "yellow" : "darkgray", x, y);
        fprintf(fp, "\\node [sChar] at (%.4f,%.4f) {\\bf \\sffamily \\textcolor{%s}{%d}};\n",
                x + 0.5, y - 0.5, i == center ? "red": "white", abs(i - center));
        
        DrawColumn(fp, context_probs[i], x, y);
    }

    // Pseudocount weights
    x = length + 1;
    y = 0.0;
    fprintf(fp, "\\filldraw [sRect,fill=yellow] (%.4f,%.4f) rectangle +(1,-1);\n", x, y);
    fprintf(fp, "\\node [sChar] at (%.4f,%.4f) {\\sffamily \\textcolor{red}{$\\boldsymbol\\nu$}};\n",
        x + 0.5, y - 0.5);
    DrawColumn(fp, &pc_probs[0], x, y);

    if (show_weights) {
        double weight_max = col_weights[0];
        for (size_t i = 1; i < length; ++i) {
            if (col_weights[i] > weight_max) {
                weight_max = col_weights[i];
            }
        }
        for (size_t i = 0; i < length; ++i) {
            fprintf(fp, "\\filldraw [sRect, draw=blue!40!black, fill=blue!40!white] "
                        "(%zu,-1) rectangle +(1,-%.4f);\n",
                    i, col_weights[i] / weight_max * weights_height);
        }
    }

    fputs("\\end{tikzpicture}\n", fp);
}

template<class Abc>
void ProbCrfStatePdfWriter<Abc>::DrawColumn(FILE* fp, double* probs, double x, double y) const {
    std::vector<size_t> res = ProbCrfStatePdfWriter<Abc>::SortRes(probs);
    for (size_t k = 0; k < res.size(); ++k) {
        size_t a = res[k];
        double height = col_height * probs[a];
        fprintf(fp, "\\filldraw [sRect, fill=%scol%c] (%.4f,%.4f) rectangle +(1,%.4f);\n",
                Abc::kName, Abc::kIntToChar[a], x, y, height);
        if (height > 1.0) {
            fprintf(fp, "\\node [sChar] at (%.4f,%.4f) {\\bf \\sffamily %c};\n",
                    x + 0.5, y + 0.5 * height, Abc::kIntToChar[a]);
        }
        y += height;
    }
}

template<class Abc>
size_t ProbCrfStatePdfWriter<Abc>::Width() const {
    return ProbCrfStatePdfWriter<Abc>::crf_state->length() + 2;
}

template<class Abc>
size_t ProbCrfStatePdfWriter<Abc>::Height() const {
    size_t height = col_height + 1;
    if (show_weights) height += weights_height;
    if (show_name) height += 1;
    return height;
}

template<class Abc>
void WeightCrfStatePdfWriter<Abc>::WriteBody(FILE* fp) {
    const CrfState<Abc>* crf_state = WeightCrfStatePdfWriter<Abc>::crf_state;
    const size_t length = crf_state->length();
    const size_t center = 0.5 * (length - 1);

    Profile<Abc> context_probs(length);
    double context_max = 0.0;
    for (size_t i = 0; i < length; ++i) {
        double pos_sum = 0.0;
        double neg_sum = 0.0;
        for (size_t a = 0; a < Abc::kSize; ++a) {
            context_probs[i][a] = exp(fabs(crf_state->context_weights[i][a]));
            if (crf_state->context_weights[i][a] > 0.0) {
                pos_sum += context_probs[i][a];
            } else {
                neg_sum += context_probs[i][a];
            }
        }
        context_max = MAX(context_max, pos_sum);
        context_max = MAX(context_max, neg_sum);
    }
    for (size_t i = 0; i < length; ++i) {
        for (size_t a = 0; a < Abc::kSize; ++a) {
            context_probs[i][a] /= context_max;
        }
    }

    ProfileColumn<Abc> pc_probs;
    double pc_max = 0.0;
    double pc_pos_sum = 0.0;
    double pc_neg_sum = 0.0;
    for (size_t a = 0; a < Abc::kSize; ++a) {
        pc_probs[a] = exp(fabs(crf_state->pc_weights[a]));
        if (crf_state->pc_weights[a] > 0.0) {
            pc_pos_sum += pc_probs[a];
        } else {
            pc_neg_sum += pc_probs[a];
        }
    }
    pc_max = MAX(pc_max, pc_pos_sum);
    pc_max = MAX(pc_max, pc_neg_sum);
    for (size_t a = 0; a < Abc::kSize; ++a) {
        pc_probs[a] /= pc_max;
    }

    fputs("\\begin{tikzpicture}[x=2.5ex, y=2.5ex]\n", fp);
    fputs("\\tikzstyle{sRect}=[rectangle,draw=black,semithick];\n", fp);
    fputs("\\tikzstyle{sChar}=[inner sep=0];\n", fp);
    fputs("\\definecolor{darkgray}{rgb}{0.3,0.3,0.3};\n", fp);

    // Name
    if (show_name) {
        fprintf(fp, "\\node [anchor=west] at (0,%.4f) {\\bf \\sffamily %s};\n",
                col_height + 1, crf_state->name.c_str());
    }

    double x;
    double y;

    // Context weights
    for (size_t i = 0; i < length; ++i) {
        y = 0.0;
        x = i;
        fprintf(fp, "\\filldraw [sRect, fill=%s] (%.4f,%.4f) rectangle +(1,-1);\n",
                i == center ? "yellow" : "darkgray", x, y);
        fprintf(fp, "\\node [sChar] at (%.4f,%.4f) {\\bf \\sffamily \\textcolor{%s}{%d}};\n",
                x + 0.5, y - 0.5, i == center ? "red": "white", abs(i - center));
        double pos_probs[Abc::kSize];
        double neg_probs[Abc::kSize];
        for (size_t a = 0; a < Abc::kSize; ++a) {
            if (crf_state->context_weights[i][a] > 0) {
                pos_probs[a] = context_probs[i][a];
                neg_probs[a] = 0.0;
            } else {
                pos_probs[a] = 0.0;
                neg_probs[a] = context_probs[i][a];
            }
        }
        DrawColumn(fp, pos_probs, x, y, false);
        DrawColumn(fp, neg_probs, x, y - 1.0, true);
    }

    // Pseudocount weights
    x = length + 1;
    y = 0.0;
    fprintf(fp, "\\filldraw [sRect,fill=yellow] (%.4f,%.4f) rectangle +(1,-1);\n", x, y);
    fprintf(fp, "\\node [sChar] at (%.4f,%.4f) {\\sffamily \\textcolor{red}{$\\boldsymbol\\nu$}};\n",
          x + 0.5, y - 0.5);
    double pos_pc_probs[Abc::kSize];
    double neg_pc_probs[Abc::kSize];
    for (size_t a = 0; a < Abc::kSize; ++a) {
        if (crf_state->pc_weights[a] > 0) {
            pos_pc_probs[a] = pc_probs[a];
            neg_pc_probs[a] = 0.0;
        } else {
            pos_pc_probs[a] = 0.0;
            neg_pc_probs[a] = pc_probs[a];
        }
    }
    DrawColumn(fp, pos_pc_probs, x, y, false);
    DrawColumn(fp, neg_pc_probs, x, y - 1, true);

    fputs("\\end{tikzpicture}\n", fp);
}

template<class Abc>
void WeightCrfStatePdfWriter<Abc>::DrawColumn(FILE* fp, double* probs, double x, double y, bool down) const {
    std::vector<size_t> res = WeightCrfStatePdfWriter<Abc>::SortRes(probs);
    for (size_t k = 0; k < res.size(); ++k) {
        size_t a = res[k];
        double height = col_height * probs[a] * (down ? -1.0 : 1.0);
        fprintf(fp, "\\filldraw [sRect, fill=%scol%c] (%.4f,%.4f) rectangle +(1,%.4f);\n",
                Abc::kName, Abc::kIntToChar[a], x, y, height);
        if (fabs(height) > 1.0) {
            fprintf(fp, "\\node [sChar] at (%.4f,%.4f) {\\bf \\sffamily %c};\n",
                    x + 0.5, y + 0.5 * height, Abc::kIntToChar[a]);
        }
        y += height;
    }
}

template<class Abc>
size_t WeightCrfStatePdfWriter<Abc>::Width() const {
    return WeightCrfStatePdfWriter<Abc>::crf_state->length() + 2;
}

template<class Abc>
size_t WeightCrfStatePdfWriter<Abc>::Height() const {
    return col_height * 2 + 1;
}

template<class Abc>
void CrfPdfWriter<Abc>::WriteBody(FILE* fp) {
    std::string s(states_per_row, 'c');
    fprintf(fp, "\\begin{tabular}{%s}\n", s.c_str());
    for (size_t i = 0; i < crf_states_.size(); ++i) {
        if (i > 0) {
            fputs(i % states_per_row == 0 ? "\\\\\n" : "&\n", fp);
        }
        if (!external_dir.empty() && !crf_states_[i].name.empty()) {
            fprintf(fp, "\\tikzsetnextfilename{%zu}\n", i + 1);
        }
        crf_state_writer_.crf_state = &crf_states_[i];
        crf_state_writer_.WriteBody(fp);
    }
    fputs("\\end{tabular}\n", fp);
}

}  // namespace cs

#endif  // CS_PDF_WRITER_INL_H_
