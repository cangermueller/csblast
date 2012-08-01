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

#include "cs.h"
#include "pdf_writer-inl.h"

namespace cs {

void PdfWriter::WriteHeader(FILE* fp) const {
    fputs("\\documentclass{article}\n", fp);
    fputs("\\usepackage[papersize={500cm,500cm},nohead,nofoot]{geometry}\n", fp);
    fputs("\\usepackage[utf8]{inputenc}\n", fp);
    fputs("\\usepackage[T1]{fontenc}\n", fp);
    fputs("\\usepackage{texshade}\n", fp);
    fputs("\\usepackage{graphicx}\n", fp);
    fputs("\\usepackage[dvipsnames]{xcolor}\n", fp);
    fputs("\\usepackage{tikz}\n", fp);
    fprintf(fp, "\\usetikzlibrary{arrows,shapes,backgrounds,petri%s}\n",
            external_dir.empty() ? "" : ",external");
    if (!external_dir.empty()) {
        fputs("\\tikzexternalize\n", fp);
        fprintf(fp, "\\tikzsetexternalprefix{%s}\n", PathCat(external_dir, ""));
    }
    fputs("\\usepackage{amsmath}\n", fp);
    fputs("\\pagestyle{empty}\n", fp);
    fputs("\\begin{document}\n", fp);
    fputs("\\definecolor{DodgerBlue}{rgb}{0.117647059,0.564705882,1.0}\n", fp);
    fputs("\\definecolor{LightSteelBlue}{HTML}{B0C4DE}\n", fp);
    fputs("\\definecolor{LightSlateGray}{HTML}{778899}\n", fp);
    fputs("\\definecolor{LightSlateBlue}{HTML}{8470FF}\n", fp);
    fputs("\\definecolor{LimeGreen}{HTML}{32CD32}\n", fp);
    fputs("\\definecolor{Lime}{HTML}{00FF00}\n", fp);
    fputs("\\definecolor{MyBlue}{rgb}{0.2549,0.4118,0.8824}\n", fp);
    fputs("\\definecolor{MyGreen}{rgb}{0.0000,0.8000,0.0000}\n", fp);

    fputs("\\definecolor{aacolW}{HTML}{00C000}\n", fp);
    fputs("\\definecolor{aacolY}{HTML}{00C000}\n", fp);
    fputs("\\definecolor{aacolF}{HTML}{00C000}\n", fp);
    fputs("\\definecolor{aacolC}{HTML}{FFFF00}\n", fp);
    fputs("\\definecolor{aacolD}{HTML}{6080FF}\n", fp);
    fputs("\\definecolor{aacolE}{HTML}{6080FF}\n", fp);
    fputs("\\definecolor{aacolL}{HTML}{02FF02}\n", fp);
    fputs("\\definecolor{aacolI}{HTML}{02FF02}\n", fp);
    fputs("\\definecolor{aacolV}{HTML}{02FF02}\n", fp);
    fputs("\\definecolor{aacolM}{HTML}{02FF02}\n", fp);
    fputs("\\definecolor{aacolK}{HTML}{FF0000}\n", fp);
    fputs("\\definecolor{aacolR}{HTML}{FF0000}\n", fp);
    fputs("\\definecolor{aacolQ}{HTML}{E080FF}\n", fp);
    fputs("\\definecolor{aacolN}{HTML}{E080FF}\n", fp);
    fputs("\\definecolor{aacolH}{HTML}{FF8000}\n", fp);
    fputs("\\definecolor{aacolP}{HTML}{A0A0A0}\n", fp);
    fputs("\\definecolor{aacolG}{HTML}{FFFFFF}\n", fp);
    fputs("\\definecolor{aacolA}{HTML}{FFFFFF}\n", fp);
    fputs("\\definecolor{aacolS}{HTML}{FFFFFF}\n", fp);
    fputs("\\definecolor{aacolT}{HTML}{FFFFFF}\n", fp);

    fputs("\\definecolor{dnacolA}{rgb}{0.0000,0.8000,0.0000}\n", fp);
    fputs("\\definecolor{dnacolC}{rgb}{0.2549,0.4118,0.8824}\n", fp);
    fputs("\\definecolor{dnacolG}{rgb}{1.0000,0.6667,0.0000}\n", fp);
    fputs("\\definecolor{dnacolT}{rgb}{0.8353,0.0000,0.0000}\n", fp);
    fputs("\\definecolor{nearlywhite}{rgb}{0.99,0.99,0.99}\n", fp);
}

void PdfWriter::WriteFooter(FILE* fp) const {
    fputs("\\end{document}\n", fp);
}

void PdfWriter::WriteToFile(std::string outfile, bool keep) {
    std::string dirname     = GetDirname(outfile);
    std::string basename    = GetBasename(outfile, false);
    std::string dirbasename = PathCat(dirname, basename);
    std::string texfile     = dirbasename + ".tex";
    std::string cmd;
    int exit_code = 0;

    // Prepare texfile
    FILE* fp = fopen(texfile.c_str(), "w");
    WriteHeader(fp);
    WriteBody(fp);
    WriteFooter(fp);
    fclose(fp);

    // Compile texfile with pdflatex
    cmd = strprintf("pdflatex %s-halt-on-error -output-directory=%s %s > /dev/null 2> /dev/null",
                    external_dir.empty() ? "" : "-shell-escape ",
                    dirname.c_str(), texfile.c_str());
    exit_code = system(cmd.c_str());
    if (exit_code) throw Exception("Error executing command '%s'!", cmd.c_str());

    // Crop resulting PDF-file with pdfcrop
    cmd = strprintf("cd %s; pdfcrop %s.pdf %s.pdf > /dev/null 2> /dev/null;",
                    dirname.c_str(), basename.c_str(), basename.c_str());
    exit_code = system(cmd.c_str());
    if (exit_code) throw Exception("Error executing command '%s'!", cmd.c_str());

    // Cleanup texfile, logfile and auxfile
    if (!keep) {
        remove(texfile.c_str());
        remove((dirbasename + ".log").c_str());
        remove((dirbasename + ".aux").c_str());
        remove((dirbasename + ".auxlock").c_str());
        if (!external_dir.empty()) {
            std::vector<std::string> files;
            GetAllFiles(external_dir, files);
            for (size_t i = 0; i < files.size(); ++i) {
                std::string ext = GetFileExt(files[i]);
                if (ext == "dpth" || ext == "log") {
                    remove(PathCat(external_dir, files[i]));
                }
            }
        }
    }
}

template<>
void AlignmentPdfWriter<AA>::WriteBody(FILE* fp) {
    FILE* fali = fopen(tempfile_.c_str(), "w");
    ali_.Write(fali, FASTA_ALIGNMENT);
    fclose(fali);

    fprintf(fp, "\\begin{texshade}{%s}\n", tempfile_.c_str());
    fputs("\\shadingmode{functional}\n", fp);
    fputs("\\funcgroup{groupWYF}{WYF}{Black}{aacolW}{upper}{bf}\n", fp);
    fputs("\\funcgroup{groupC}{C}{Black}{aacolC}{upper}{bf}\n", fp);
    fputs("\\funcgroup{groupDE}{DE}{Black}{aacolD}{upper}{bf}\n", fp);
    fputs("\\funcgroup{groupLIVM}{LIVM}{Black}{aacolL}{upper}{bf}\n", fp);
    fputs("\\funcgroup{groupKR}{KR}{Black}{aacolK}{upper}{bf}\n", fp);
    fputs("\\funcgroup{groupQN}{QN}{Black}{aacolQ}{upper}{bf}\n", fp);
    fputs("\\funcgroup{groupH}{H}{Black}{aacolH}{upper}{bf}\n", fp);
    fputs("\\funcgroup{groupP}{P}{Black}{aacolP}{upper}{bf}\n", fp);
    fputs("\\funcgroup{groupGAST}{GAST}{Black}{lightgray}{upper}{bf}\n", fp);
    // fputs("\\funcgroup{groupGAST}{GAST}{Black}{white}{upper}{bf}\n", fp);
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
    // fputs("\\hidenames\n", fp);
    fputs("\\end{texshade}\n", fp);
}

}  // namespace cs
