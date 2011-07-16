// Copyright 2010, Andreas Biegert

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
    fputs("\\usetikzlibrary{arrows,shapes,backgrounds,petri}\n", fp);
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

void PdfWriter::WriteToFile(std::string outfile) {
    std::string dirname  = GetDirname(outfile);
    std::string basename = GetBasename(outfile, false);
    std::string texfile  = dirname + kDirSep + basename + ".tex";
    std::string logfile  = dirname + kDirSep + basename + ".log";
    std::string auxfile  = dirname + kDirSep + basename + ".aux";
    std::string cmd;
    int exit_code = 0;

    // Prepare texfile
    FILE* fp = fopen(texfile.c_str(), "w");
    WriteHeader(fp);
    WriteBody(fp);
    WriteFooter(fp);
    fclose(fp);

    // Compile texfile with pdflatex
    //  cmd = strprintf("pdflatex -output-directory=%s %s > /dev/null 2> /dev/null",
    //                dirname.c_str(), texfile.c_str());
    cmd = strprintf("pdflatex -output-directory=%s %s > /dev/null 2> /dev/null",
                    dirname.c_str(), texfile.c_str());
    exit_code = system(cmd.c_str());
    if (exit_code) throw Exception("Error executing command '%s'!", cmd.c_str());

    // Crop resulting PDF-file with pdfcrop
    cmd = strprintf("cd %s; pdfcrop %s.pdf %s.pdf > /dev/null 2> /dev/null;",
                    dirname.c_str(), basename.c_str(), basename.c_str());
    exit_code = system(cmd.c_str());
    if (exit_code) throw Exception("Error executing command '%s'!", cmd.c_str());

    // Cleanup texfile, logfile and auxfile
    remove(texfile.c_str());
    remove(logfile.c_str());
    remove(auxfile.c_str());
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
