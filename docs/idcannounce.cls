%% $Id: idcannounce.cls 121 2012-10-11 09:36:46Z aschmidt $
%% Copyright 2009 Armin C. Schmidt <dunkelschorsch AT googlemail DOT com>
%
% This work may be distributed and/or modified under the
% conditions of the LaTeX Project Public License, either version 1.3
% of this license or (at your option) any later version.
% The latest version of this license is in
%   http://www.latex-project.org/lppl.txt
% and version 1.3 or later is part of all distributions of LaTeX
% version 2005/12/01 or later.
%
% This work has the LPPL maintenance status `maintained'.
%
% The Current Maintainer of this work is Armin C. Schmidt.
%
% This work consists of all files listed in manifest.txt
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{idcannounce}[2009/04/01 LaTeX2e class for announcements at the IDC]

\RequirePackage{ifxetex}
\ifxetex
  \RequirePackage{fontspec}
  \setsansfont[ExternalLocation = fonts/,
    Extension = .ttf,
    UprightFont = * Roman,
    BoldFont = * Bold,
    ItalicFont = * Italic,
    BoldItalicFont = * Bold Italic]{Humnst777 BT}
    \let\familydefault\sfdefault
\else
  \RequirePackage[T1]{fontenc}
  \RequirePackage{lmodern}
  \AtBeginDocument{\usepackage{microtype}}
\fi

\RequirePackage{tikz}
\RequirePackage{tabularx}

\definecolor{lmkblue}{RGB}{47,59,154}
\definecolor{wavegreen}{cmyk}{1.0, 0.0, 1.0, 0.0}

\newlength{\idc@hmargin}
\setlength{\idc@hmargin}{22mm}

\newlength{\idc@topmargin}
\setlength{\idc@topmargin}{20mm}

\newlength{\idc@logo@height}
\setlength{\idc@logo@height}{27mm}
\pgfdeclareimage[height=\idc@logo@height]{idc@logo}{graphics/idc_header}

\newcommand{\idc@head}{
  \thispagestyle{empty}
  \begin{tikzpicture}[overlay,remember picture]
    \tikzstyle{every node}=[text centered]
% Logo
    \path (current page.north)+(0,-6mm) node[anchor=north] (logo) {\pgfuseimage{idc@logo}};
  \end{tikzpicture}
  \vspace*{0.3\idc@logo@height}\par
}

\AtBeginDocument{\idc@head}

\LoadClass[12pt,a4paper]{scrartcl}

\setlength{\parindent}{0mm}
\setlength{\parskip}{3mm}

\RequirePackage[paper=a4paper,hmargin=\idc@hmargin,top=\idc@topmargin,vmarginratio=1:1]{geometry}
\addtokomafont{section}{\large}

\RequirePackage{titlesec}
\titlespacing{\section}{0em}{0.5em}{0em}

\newcommand{\signature}[2]{%
 \begin{tabularx}{\textwidth}{@{}lX>{\centering\let\newline\\\arraybackslash\hspace{0pt}}p{65mm}}
  Supervisor: #1 & ~ & \dotfill \\[-1mm]
  ~ & ~ & \small (#2)
 \end{tabularx}
}

\newcommand*{\heading}[2][]{
  {\centering
    \usekomafont{disposition}%
    \if#1\@empty%
      \LARGE\textcolor{wavegreen}{#2}%
    \else%
      \LARGE\textcolor{wavegreen}{#1}\\[1ex]%
      \large#2%
    \fi\\[1.5\baselineskip]}
}
