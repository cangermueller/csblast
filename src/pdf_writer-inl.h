// Copyright 2010, Andreas Biegert

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
POHmmPdfWriter<Abc>::POHmmPdfWriter(const POHmm<Abc>& hmm,
                                    const ContextLibrary<Abc>* lib,
                                    const IndexPath* path)
        : hide_numbering(false),
          hide_profiles(false),
          hide_alignment(false),
          hide_paths(false),
          reverse(false),
          pmin(0.02),
          hmm_(hmm),
          lib_(lib),
          path_(path),
          grid_(hmm) {
    // Generate temporary filename as basename for alignment block PDFs
    char tempfilename[] = "/tmp/alignment_XXXXXX";
    const int captured_fd = mkstemp(tempfilename);
    if (!captured_fd) throw Exception("Unable to create unique alignment filename!");
    tempfile_ = tempfilename;
}

template<class Abc>
POHmmPdfWriter<Abc>::~POHmmPdfWriter() {
    remove(tempfile_.c_str());
    for (size_t b = 1; b < grid_.blocks.size() - 1; ++b)
        remove(alignment_filename(b).c_str());
}

template<class Abc>
void POHmmPdfWriter<Abc>::WriteBody(FILE* fp) {
    const char kStartEndColor[] = "LightSteelBlue";
    const Graph& g = hmm_.g;
    const Vertex kEndVertex = hmm_.vertices.size();  // virtual END vertex
    double s = reverse ? -1.0 : 1.0;  // draw PO-HMM in reverse from right to left
    double z = 0.0;  // y-coordinate baseline

    DefineContextColors(fp);
    fputs("\\sffamily\n", fp);
    fprintf(fp, "\\begin{tikzpicture}[x=%.5fcm,y=%.5fcm,\n", kResWidth, kScaleY);
    fputs("vertex/.style={circle,thin,inner sep=0pt,minimum size=1.0mm},\n", fp);
    fputs("skipedge/.style={to path={-- ++(0,#1) -| (\\tikztotarget)}},\n", fp);
    fputs("logo/.style={rectangle,draw=black,inner sep=0pt,minimum size=0mm,"
          "anchor=south,semithick},\n", fp);
    fputs("alibox/.style={rectangle,thin,inner sep=0pt,draw=black}]\n", fp);

    // Draw edges between vertices first
    EdgeIter ei, ei_end;
    double hmax = 1, hmin = -1;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        double sat = 20 + 80 * g[*ei].weight;
        Vertex u = source(*ei, g);
        Vertex v = target(*ei, g);
        if (v == kStartEndVertex) v = kEndVertex;

        if ((v != kEndVertex && hmm_.indices[u] == hmm_.indices[v] - 1) ||
            (v == kEndVertex && hmm_.indices[u] == hmm_.size())) {
            // Draw flat edge from left to right
            fprintf(fp, "\\path (%.2f,0) edge[semithick,draw=black!%.0f] "
                    "(%.2f,0);\n", s * grid_.coords[u], sat, s * grid_.coords[v]);
        } else {
            // Draw skip-edge from left to right with height dependant on distance
            double h = pow(fabs(grid_.coords[u] - grid_.coords[v]), 0.7);
            if ((g[u].alis & g[v].alis) == g[u].alis) h *= -1;
            hmax = MAX(h, hmax);
            hmin = MIN(h, hmin);
            double x = s * grid_.coords[u] + s * 0.1;
            double y = s * grid_.coords[v] - s * 0.1;
            fprintf(fp, "\\path (%.2f,0) edge[semithick,skipedge=%.2f,line join=round,"
                    "draw=black!%.0f] (%.2f,0);\n", x, h, sat, y);
        }
    }

    // Draw vertex circles
    fprintf(fp, "\\node[vertex,draw=black,fill=%s] at (%.2f,%.2f) (v%zu) {};\n",
            kStartEndColor, s * grid_.coords[kStartEndVertex], z, kStartEndVertex);
    for (ConstVertexVecIter v = hmm_.begin() + 1; v != hmm_.end(); ++v) {
        int fill = g[*v].neff < hmm_.neff ? 40 : 100;
        fprintf(fp, "\\node[vertex,draw=black,fill=black!%d] at (%.2f,%.2f) "
                "(v%zu) {};\n", fill, s * grid_.coords[*v], z, *v);
    }
    fprintf(fp, "\\node[vertex,draw=black,fill=%s] at (%.2f,%.2f) (v%zu) {};\n",
            kStartEndColor, s * grid_.coords[kEndVertex], z, kEndVertex);
    z += hmax;

    if (!hide_numbering) {
        // Draw vertex numbers
        double rotate = s * 90;
        std::string anchor = reverse ? "east" : "west";
        for (ConstVertexVecIter v = hmm_.begin() + 1; v != hmm_.end(); ++v)
            fprintf(fp, "\\node[rotate=%.0f,anchor=%s] at (%.2f,%.2f) {\\tiny %zu};\n",
                    rotate, anchor.c_str(), s * grid_.coords[*v], z, *v);
        z += 2.0 + 1.6 * strprintf("%zu", hmm_.size()).length();
    }

    if (!hide_alignment) {
        // Draw alignment blocks
        GenerateAlignmentBlocks();
        for (size_t b = 1; b < grid_.blocks.size() - 1; ++b) {
            double x = s * grid_.coords[grid_.blocks[b].front()] - 1.15 * s;
            std::string anchor = reverse ? "south east" : "south west";
            fprintf(fp, "\\node[anchor=%s] at (%.2f,%.2f) {\\includegraphics{%s}};\n",
                    anchor.c_str(), x, z, alignment_filename(b).c_str());
        }
        z += 2.5 + 5.0 * hmm_.seqs.size();
    }

    if (lib_ != NULL && !hide_profiles) {
        // Draw context-state profiles
        for (ConstVertexVecIter v = hmm_.begin() + 1; v != hmm_.end(); ++v) {
            // Sort contexts by color and probability
            std::vector<ColorProbPair> pairs;
            for (ContextIter c = g[*v].contexts.begin(); c != g[*v].contexts.end(); ++c)
                if (c->prob > pmin)
                    pairs.push_back(std::make_pair((*lib_)[c->index].color, c->prob));
            std::sort(pairs.begin(), pairs.end());
            std::reverse(pairs.begin(), pairs.end());
            // Draw profile logos
            double x = s * grid_.coords[*v] - 0.5 * s;
            double y = z;
            for (size_t k = 0; k < pairs.size(); ++k) {
                double red   = 100 * pairs[k].first.red;
                double green = 100 * pairs[k].first.green;
                double blue  = 100 * pairs[k].first.blue;
                double prob  = pairs[k].second;
                fprintf(fp, "\\filldraw [logo,fill=col%.0f%.0f%.0f] (%.4f,%.4f) "
                        "rectangle +(%.4f,%.4f);\n", red, green, blue, x, y, s, prob * 20);
                y += prob * 20;
            }
        }
        z += 20;
    } else if (!hide_profiles) {
        // Draw alphabet profiles
        for (ConstVertexVecIter v = hmm_.begin() + 1; v != hmm_.end(); ++v) {
            // Sort residues by functional group and probability
            std::vector<GroupProbResTuple> tuples;
            for (size_t a = 0; a < Abc::kSize; ++a)
                if (g[*v].probs[a] > pmin)
                    tuples.push_back(boost::make_tuple(-Abc::kFuncGroup[a], g[*v].probs[a], a));
            std::sort(tuples.begin(), tuples.end());
            std::reverse(tuples.begin(), tuples.end());
            // Draw profile logos
            double x = s * grid_.coords[*v] - 0.5 * s;
            double y = z;
            for (size_t k = 0; k < tuples.size(); ++k) {
                double prob  = boost::get<1>(tuples[k]);
                uint8_t res  = boost::get<2>(tuples[k]);
                fprintf(fp, "\\filldraw[logo,fill=%scol%c] (%.4f,%.4f) rectangle +(%.4f,%.4f);\n",
                        Abc::kName, Abc::kIntToChar[res], x, y, s, prob * 20);
                if (prob > 0.1)
                    fprintf(fp, "\\node at (%.4f,%.4f) {\\tiny \\bf \\sffamily %c};\n",
                            x + s * 0.5, y + prob * 10, Abc::kIntToChar[res]);
                y += prob * 20;
            }
        }
        z += 20;
    }

    if (!hide_paths) {
        // Draw alignment paths
        for (size_t n = 0; n < hmm_.num_alis; ++n) {
            Vertex u = kStartEndVertex;
            bool on = false;
            for (ConstVertexVecIter vi = hmm_.begin() + 1; vi != hmm_.end(); ++vi) {
                if (!on && g[*vi].alis.test(n)) {
                    // Found start of alignment path => remember start vertex
                    u = *vi;
                    on = true;
                    if (vi + 1 != hmm_.end()) continue;
                }
                if (on && (vi + 1 == hmm_.end() || !g[*vi].alis.test(n))) {
                    // Found end of alignment path => print color box
                    Vertex v = vi + 1 == hmm_.end() ? *vi : *(vi - 1);
                    double x1 = s * grid_.coords[u] - 0.5 * s;
                    double x2 = s * grid_.coords[v] + 0.5 * s;
                    double y = hmin - n - 2;
                    fprintf(fp, "\\filldraw[alibox,fill=%s] (%.2f,%.2f) rectangle "
                            "(%.2f,%.2f);\n", kAliColors[n], x1, y, x2, y - 1);
                    on = false;
                }
            }
        }
    }

    fputs("\\begin{pgfonlayer}{background}\n", fp);
    // For each block draw yellow boxes as background
    for (size_t b = 0; b < grid_.blocks.size(); ++b) {
        Vertex u = reverse ? grid_.blocks[b].back() : grid_.blocks[b].front();
        Vertex v = reverse ? grid_.blocks[b].front() : grid_.blocks[b].back();
        double flow = u == kEndVertex ? 1.0 : g[u].flow;
        fprintf(fp, "\\filldraw [line width=1mm,Goldenrod!%.0f,join=round] "
                "(v%zu.north -| v%zu.west) rectangle (v%zu.south -| v%zu.east);\n",
                20 + 80 * flow, u, u, v, v);
    }
    // Highlight vertex path if it was given
    if (path_ != NULL) {
        std::string colorMM = PairStateToColor(MM);
        fprintf(fp, "\\node[vertex,draw=%s,fill=%s,,minimum size=1.5mm] "
                "at (0,0) {};\n", colorMM.c_str(), colorMM.c_str());
        for (IndexPath::const_iterator p = path_->begin(); p != path_->end(); ++p) {
            std::string color = PairStateToColor(p->second);
            fprintf(fp, "\\node[vertex,draw=%s,fill=%s,,minimum size=1.5mm] "
                    "at (%.2f,0) {};\n", color.c_str(), color.c_str(),
                    s * grid_.coords[p->first]);
        }
        fprintf(fp, "\\node[vertex,draw=%s,fill=%s,,minimum size=1.5mm] "
                "at (%.2f,0) {};\n", colorMM.c_str(), colorMM.c_str(),
                s * grid_.coords[kEndVertex]);
    }
    fputs("\\end{pgfonlayer}\n", fp);
    fputs("\\end{tikzpicture}\n", fp);
}

template<class Abc>
void POHmmPdfWriter<Abc>::GenerateAlignmentBlocks() const {
    const Graph& g = hmm_.g;

    for (size_t b = 1; b < grid_.blocks.size() - 1; ++b) {
        Alignment<Dna> ali(grid_.blocks[b].size(), hmm_.seqs.size());
        for (size_t i = 0; i < grid_.blocks[b].size(); ++i) {
            Vertex v = grid_.blocks[b][i];
            size_t j = reverse ? grid_.blocks[b].size() - 1 - i : i;
            for (ConstColIter r = g[v].col.begin(); r != g[v].col.end(); ++r)
                ali[j][r->first] = hmm_.seqs[r->first][r->second];
        }
        AlignmentPdfWriter<Abc> ali_writer(ali);
        ali_writer.WriteToFile(alignment_filename(b));
    }
}

template<class Abc>
void POHmmPdfWriter<Abc>::DefineContextColors(FILE* fp) const {
    if (lib_ == NULL) return;
    std::set<Color> colors;
    for (size_t k = 0; k < lib_->size(); ++k)
        colors.insert((*lib_)[k].color);
    for (std::set<Color>::iterator c = colors.begin(); c != colors.end(); ++c)
        fprintf(fp, "\\definecolor{col%.0f%.0f%.0f}{rgb}{%.2f,%.2f,%.2f}\n",
                100 * c->red, 100 * c->green, 100 * c->blue, c->red, c->green,
                c->blue);
}

// For converting from ASCII to the amino acid code
template<class Abc>
const char* POHmmPdfWriter<Abc>::kAliColors[] = {
    "DodgerBlue", "Dandelion", "OrangeRed", "LimeGreen", "Periwinkle", "Apricot",
    "SeaGreen", "RedOrange",  "RoyalBlue", "Gray", "Bittersweet", "Magenta"
};


template<class Abc>
PosteriorMatrixPdfWriter<Abc>::PosteriorMatrixPdfWriter(
        const AlignmentMatrices<Abc>& data,
        const ContextLibrary<Abc>* lib,
        const PairAlignment* ali)
        : data_(data),
          lib_(lib),
          ali_(ali),
          grid_x_(data.x),
          grid_y_(data.y) {
    // Generate temporary filename as basename for alignment block PDFs
    char tempfilename[] = "/tmp/alignment_XXXXXX";
    const int captured_fd = mkstemp(tempfilename);
    if (!captured_fd) throw Exception("Unable to create unique alignment filename!");
    tempfile_ = tempfilename;
}

template<class Abc>
PosteriorMatrixPdfWriter<Abc>::~PosteriorMatrixPdfWriter() {
    remove(tempfile_.c_str());
    remove((tempfile_ + "_x.pdf").c_str());
    remove((tempfile_ + "_y.pdf").c_str());
}

template<class Abc>
void PosteriorMatrixPdfWriter<Abc>::WriteBody(FILE* fp) {
    const Vertex kEndX = data_.x.vertices.size();  // virtual END vertex in x
    const Vertex kEndY = data_.y.vertices.size();  // virtual END vertex in y

    // Generate PO-HMM PDFs for x and y axis
    scoped_ptr<IndexPath> path_x;
    if (ali_ != NULL) path_x.reset(new IndexPath(GetPathInX(*ali_)));
    POHmmPdfWriter<Abc> x_writer(data_.x, lib_, path_x.get());
    std::string hmmfile_x(tempfile_ + "_x.pdf");
    x_writer.reverse = true;
    x_writer.hide_profiles = true;
    x_writer.hide_alignment = false;
    x_writer.hide_paths = false;
    x_writer.WriteToFile(hmmfile_x);

    scoped_ptr<IndexPath> path_y;
    if (ali_ != NULL) path_y.reset(new IndexPath(GetPathInY(*ali_)));
    POHmmPdfWriter<Abc> y_writer(data_.y, lib_, path_y.get());
    std::string hmmfile_y(tempfile_ + "_y.pdf");
    y_writer.hide_profiles = true;
    y_writer.hide_alignment = false;
    y_writer.hide_paths = false;
    y_writer.WriteToFile(hmmfile_y);

    // Write TikZ picture
    fprintf(fp, "\\begin{tikzpicture}[x=%.5fcm,y=%.5fcm,\n", kResWidth, kResWidth);
    fprintf(fp, "cell/.style={rectangle,very thin,draw=lightgray,inner sep=0pt, "
            "minimum size=%.5fcm},\n", kResWidth);
    fprintf(fp, "path/.style={rectangle,semithick,inner sep=0pt, "
            "minimum size=%.5fcm}]\n", kResWidth - 0.022);

    // Add axis labels
    fprintf(fp, "\\node[anchor=south east,inner sep=0mm, line width=0mm, rotate=90] "
            "at (-1.0,0.60) {\\includegraphics{%s}};\n", hmmfile_x.c_str());
    fprintf(fp, "\\node[anchor=south west,inner sep=0mm, line width=0mm] "
            "at (-0.65,1.0) {\\includegraphics{%s}};\n", hmmfile_y.c_str());

    // Draw posterior matrix cells shaded by probability
    std::string pcol;
    for (ConstVertexVecIter xi = data_.x.begin(); xi != data_.x.end(); ++xi) {
        for (ConstVertexVecIter yj = data_.y.begin(); yj != data_.y.end(); ++yj) {
            size_t i = *xi;
            size_t j = *yj;
            double x =  grid_y_.coords[j];
            double y = -grid_x_.coords[i];
            int red   = iround(100.0 * data_.P[i][j].MM);
            int green = iround(100.0 * data_.P[i][j].MI);
            int blue  = iround(100.0 * data_.P[i][j].IM);
            int white = 100 - red - green - blue;
            pcol = strprintf("{rgb:red,%d;green,%d;blue,%d;white,%d}",
                             red, green, blue, white);
            fprintf(fp, "\\node[cell,fill=%s] at (%.2f,%.2f) {};\n", pcol.c_str(), x, y);
        }
    }
    // Draw END-END cell
    fprintf(fp, "\\node[cell,fill=red] at (%d,%d) {};\n",
            grid_y_.coords[kEndY], -grid_x_.coords[kEndX]);

    // Highlight alignment path if a pairwise alignment was given
    std::string pathcol = "black";
    if (ali_ != NULL) {
        size_t i_last = kStartEndVertex, j_last = kStartEndVertex;
        // Draw START-START cell
        fprintf(fp, "\\node[path,draw=%s] at (0,0) (i0j0) {};\n", pathcol.c_str());
        for (ConstPathIter s = ali_->path.begin(); s != ali_->path.end(); ++s) {
            double x =  grid_y_.coords[s->j];
            double y = -grid_x_.coords[s->i];
            fprintf(fp, "\\node[path,draw=%s] at (%.2f,%.2f) (i%zuj%zu) {};\n",
                    pathcol.c_str(), x, y, s->i, s->j);
            if ((grid_x_.coords[s->i] - grid_x_.coords[i_last] > 1) ||
                (grid_y_.coords[s->j] - grid_y_.coords[j_last] > 1)) {
                fprintf(fp, "\\draw[->,draw=%s,semithick] (i%zuj%zu) to node {} "
                        "(i%zuj%zu);\n", pathcol.c_str(), s->i, s->j, i_last, j_last);
            }
            i_last = s->i;
            j_last = s->j;
        }
        fprintf(fp, "\\node[path,draw=%s] at (%d,%d) (i%zuj%zu) {};\n", pathcol.c_str(),
                grid_y_.coords[kEndY], -grid_x_.coords[kEndX], kEndX, kEndY);
        fprintf(fp, "\\draw[->,draw=%s,thick] (i%zuj%zu) to node {} (i%zuj%zu);\n",
                pathcol.c_str(), kEndX, kEndY, i_last, j_last);
    }

    fputs("\\end{tikzpicture}\n", fp);
}


template<class AS, class Abc>
AbstractStateMatrixPdfWriter<AS, Abc>::AbstractStateMatrixPdfWriter(
        const AbstractStateMatrix<AS>& matrix,
        const ContextLibrary<Abc>& contexts,
        const ContextLibrary<Abc>& alphabet)
        : score_min(-8.0),
          score_max(8.0),
          matrix_(matrix),
          contexts_(contexts),
          alphabet_(alphabet) {
    assert_eq(matrix_.num_contexts(), contexts_.size());
    // Generate temporary filename as basename for alignment block PDFs
    char tempfilename[] = "/tmp/abstract_state_alphabet_XXXXXX";
    const int captured_fd = mkstemp(tempfilename);
    if (!captured_fd) throw Exception("Unable to create unique alignment filename!");
    tempfile_ = tempfilename;
}

template<class AS, class Abc>
AbstractStateMatrixPdfWriter<AS, Abc>::~AbstractStateMatrixPdfWriter() {
    remove(tempfile_.c_str());
    remove((tempfile_ + ".pdf").c_str());
}

template<class AS, class Abc>
Color AbstractStateMatrixPdfWriter<AS, Abc>::GetColor(size_t k, size_t a) const {
    Color col;
    if (matrix_.q(k,a) == 0.0) return col;
    if (matrix_.s(k,a) >= 0.0) {
        double s = MIN(1.0, matrix_.s(k,a) / score_max);
        col.red   = s;
        col.green = 1.0 - s;
        col.blue  = 0.0;
    } else {
        double s = MIN(1.0, -matrix_.s(k,a) / -score_min);
        col.red   = 0.0;
        col.green = 1.0 - s;
        col.blue  = s;
    }
    return col;
}

template<class AS, class Abc>
void AbstractStateMatrixPdfWriter<AS, Abc>::DefineContextColors(FILE* fp) const {
    std::set<Color> colors;
    for (size_t k = 0; k < contexts_.size(); ++k)
        colors.insert(contexts_[k].color);
    for (std::set<Color>::iterator c = colors.begin(); c != colors.end(); ++c)
        fprintf(fp, "\\definecolor{col%.0f%.0f%.0f}{rgb}{%.2f,%.2f,%.2f}\n",
                100 * c->red, 100 * c->green, 100 * c->blue,
                c->red, c->green, c->blue);
}

template<class AS, class Abc>
void AbstractStateMatrixPdfWriter<AS, Abc>::WriteBody(FILE* fp) {
    // Generate PDF of abstract state alphabet
    ContextLibraryPdfWriter<Abc> alphabet_writer(alphabet_);
    alphabet_writer.margin = 0;
    alphabet_writer.ncols  = alphabet_.size();
    std::string alphabet_pdf = tempfile_ + ".pdf";
    alphabet_writer.WriteToFile(alphabet_pdf);

    // Write TikZ picture
    DefineContextColors(fp);
    fprintf(fp, "\\begin{tikzpicture}[x=%.5fcm,y=%.5fcm,\n", kResWidth, kResWidth);
    fprintf(fp, "cname/.style={rectangle,thin,draw=black,inner sep=0.2mm,minimum "
            "size=%.5fcm,anchor=east,rounded corners=1pt},\n", kResWidth);
    fprintf(fp, "cell/.style={rectangle,very thin,draw=white,inner sep=0pt, "
            "minimum size=%.5fcm}]\n", kResWidth);

    // Show abstract state alphabet above score matrix
    fprintf(fp, "\\node[anchor=south west,inner sep=0mm, line width=0mm] "
            "at (-0.65,1.0) {\\includegraphics{%s}};\n", alphabet_pdf.c_str());

    // Draw substition matrix cells colored by scores
    const size_t kmax = MIN(1500, static_cast<int>(contexts_.size()));
    std::string str;
    for (size_t k = 0; k < kmax; ++k) {
        // Draw name badge of context k
        const ContextProfile<Abc>& cp = contexts_[k];
        fprintf(fp, "\\message{%zu %s}\n", k, cp.name.c_str());
        fprintf(fp, "\\node[cname,fill=col%.0f%.0f%.0f] at (-1,-%zu) "
                "{\\tiny \\bf \\sffamily %s};\n", 100 * cp.color.red,
                100 * cp.color.green, 100 * cp.color.blue, k,
                cp.name.empty() ? strprintf("%zu", k).c_str() : cp.name.c_str());

        for (size_t a = 0; a < AS::kSize; ++a) {
            Color col = GetColor(k, a);
            int red   = iround(100.0 * col.red);
            int green = iround(100.0 * col.green);
            int blue  = iround(100.0 * col.blue);
            str = strprintf("{rgb:red,%d;green,%d;blue,%d}", red, green, blue);
            fprintf(fp, "\\node[cell,fill=%s] at (%zu,-%zu) {};\n", str.c_str(), a, k);
        }
    }

    fputs("\\end{tikzpicture}\n", fp);
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
        crf_state_writer_.crf_state = &crf_states_[i];
        crf_state_writer_.WriteBody(fp);
        fputs((i + 1) == states_per_row ? "\\\\" : "&", fp);
    }
    fputs("\\end{tabular}\n", fp);

    /*
    double x = 0.0;
    double y = 0.0;
    double sep = 1;
    fputs("\\begin{tikzpicture}[x=2.5ex,y=2.5ex]\n", fp);
    for (size_t i = 0; i < crf_states_.size(); ++i) {
        fprintf(fp, "\\node at (%.4f,%.4f) {\n", x, y);
        crf_state_writer_.crf_state = &crf_states_[i];
        crf_state_writer_.WriteBody(fp);
        fputs("};\n", fp);
        if ((i + 1) % states_per_row == 0) {
            x = 0.0;
            y -= crf_state_writer_.Height() + sep;
        } else {
            x += crf_state_writer_.Width() + sep;
        }
    }
    fputs("\\end{tikzpicture}\n", fp);
    */
}

}  // namespace cs

#endif  // CS_PDF_WRITER_INL_H_
