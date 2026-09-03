\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{reviewresponse}[2025/01/16 v2.0.1 class for writing rebuttal letters for reviews of submitted papers]
% Author: Karl-Ludwig Besser
% Email: karl-ludwig.besser@liu.se

\RequirePackage{kvoptions}

\newif\ifreviewresponselevel
\reviewresponseleveltrue

\newif\ifreviewresponsecompleted

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{scrartcl}}
\ProcessOptions\relax
\LoadClass{scrartcl}


\RequirePackage[breakable,skins]{tcolorbox}
\RequirePackage{xcolor}
\RequirePackage{xparse}


\newcommand*{\journal}[1]{\renewcommand*{\@journal}{#1}}
\newcommand*{\@journal}{}
\newcommand*{\editorname}[1]{\renewcommand*{\@editorname}{#1}}
\newcommand*{\@editorname}{}
\newcommand*{\manuscript}[1]{\renewcommand*{\@manuscript}{#1}}
\newcommand*{\@manuscript}{}

\newcommand*{\thejournal}{\@journal}
\newcommand*{\thetitle}{\@title}
\newcommand*{\theauthor}{\@author}
\newcommand*{\thedate}{\@date}
\newcommand*{\theeditor}{\@editorname}
\newcommand*{\themanuscript}{\@manuscript}

\newcommand{\listprioH}{}
\newcommand{\listprioM}{}
\newcommand{\listprioL}{}

%%% Colors
% \definecolor{colorcommentfg}{RGB}{0,63,87}
\definecolor{colorcommentfg}{RGB}{0,0,0}
% \definecolor{colorcommentbg}{HTML}{e0f0f6}
\definecolor{colorcommentbg}{HTML}{f0f0f0}
% \definecolor{colorcommentframe}{RGB}{0,112,155}
\definecolor{colorcommentframe}{RGB}{0,0,0}

\colorlet{colorchangebg}{black!2}
\colorlet{colorchangeframe}{black!20}

\colorlet{colorconcludingresponsefg}{black!80}
\colorlet{colorconcludingresponsebg}{black!2}
\colorlet{colorconcludingresponseframe}{black!20}
%%%


%%% General Settings
\setlength{\parindent}{0pt}
\setlength{\parskip}{.3em plus .2em minus .3em}
\renewcommand{\baselinestretch}{1.5}

\widowpenalty10000
\clubpenalty10000
%%%

%%% Counters
\newcounter{editorresponse}
\setcounter{editorresponse}{0}
\newcounter{reviewer}
\setcounter{reviewer}{0}
\newcounter{revcomment}[reviewer]
\setcounter{revcomment}{0}
%%%

%%% Commands
\renewcommand*{\maketitle}{%
	\begin{titlepage}
		\begin{center}
			\vspace*{1cm}
			\large{Responses to Reviewers' Comments for Manuscript \@manuscript}\\
			
			\vspace{1cm}
			\textbf{\LARGE{\@title}}\\
			\vspace{1cm}
			\large{Addressed Comments for Publication to}
			
			\Large{\@journal}
			
			\large{by}
			
			{\large{\@author}}
		\end{center}
	\end{titlepage}
}

\newcommand*{\generalhowto}{
	\clearpage
	% \@ifundefined{pdfbookmark}{}{%
	% 	\pdfbookmark[1]{Reponse to the Editor}{hyperref@editor}%
	% }%
	\section*{General Approach of How to Answering Comments}
}

\newcommand*{\reviewresponse@editorlabel}{Editor}
\newcommand*{\reviewresponse@reviewerlabel}{Reviewer~\arabic{reviewer}}

\newcommand*{\editor}{
	\clearpage
	\refstepcounter{editorresponse}%
	\addcontentsline{toc}{section}{Author's Response to \reviewresponse@editorlabel} % Add to TOC
	\@ifundefined{pdfbookmark}{}{%
		\pdfbookmark[1]{Response to \reviewresponse@editorlabel}{hyperref@editor\arabic{editorresponse}}%
	}%
	\section*{Authors' Response to \reviewresponse@editorlabel}
}

\newcommand*{\reviewer}{
	\clearpage
	\refstepcounter{reviewer}%
	\addcontentsline{toc}{section}{\reviewresponse@reviewerlabel} % Add to TOC
	\@ifundefined{pdfbookmark}{}{%
		\pdfbookmark[1]{\reviewresponse@reviewerlabel}{hyperref@reviewer\arabic{reviewer}}%
	}%
	\section*{Authors' Response to \reviewresponse@reviewerlabel}
}
%%%

%%% Blocks %%%
\newenvironment{generalcomment}[1][General Comments]{%
	% \addcontentsline{toc}{subsection}{#1}%
	\begin{tcolorbox}[attach title to upper, title={#1}, after title={.\enskip}, fonttitle={\rmfamily\bfseries}, coltitle={colorcommentfg}, colback={colorcommentbg}, colframe={colorcommentframe},bookmark={#1}]
	}{\end{tcolorbox}}


% colors for priority levels
\definecolor{prio1}{HTML}{D32F2F} % High: Red
\definecolor{prio2}{HTML}{F57C00} % Medium: Orange
\definecolor{prio3}{HTML}{1976D2} % Low: Blue

\newcommand*{\reviewresponse@prioritytitle}{ }
\newcommand*{\reviewresponse@priorityframe}{colorchangeframe}
\newcommand*{\reviewresponse@prioritybg}{colorchangebg}

\newcommand*{\reviewresponse@setprioritystyle}[1]{%
	\ifreviewresponselevel
		\ifnum#1=1
			\renewcommand*{\reviewresponse@prioritytitle}{High Priority}%
			\colorlet{\reviewresponse@priorityframe}{prio1}%
			\colorlet{\reviewresponse@prioritybg}{prio1!5}%
		\else\ifnum#1=2
			\renewcommand*{\reviewresponse@prioritytitle}{Medium Priority}%
			\colorlet{\reviewresponse@priorityframe}{prio2}%
			\colorlet{\reviewresponse@prioritybg}{prio2!5}%
		\else
			\renewcommand*{\reviewresponse@prioritytitle}{Low Priority}%
			\colorlet{\reviewresponse@priorityframe}{prio3}%
			\colorlet{\reviewresponse@prioritybg}{prio3!5}%
		\fi\fi
	\else
		\renewcommand*{\reviewresponse@prioritytitle}{}
		\colorlet{\reviewresponse@priorityframe}{colorchangeframe}
		\colorlet{\reviewresponse@prioritybg}{colorchangebg}
	\fi
}

\newcommand*{\reviewresponse@setcompletedstyle}[1]{%
	\reviewresponsecompletedfalse
	\ifnum#1=1
		\reviewresponsecompletedtrue
	\fi
	\ifreviewresponsecompleted
		\colorlet{\reviewresponse@priorityframe}{black!30}
		\colorlet{\reviewresponse@prioritybg}{black!6}
	\fi
}

\renewcommand{\therevcomment}{\the\value{reviewer}.\the\value{revcomment}}

\NewDocumentEnvironment{revcomment}{ O{3} O{0} m }{% 
    \refstepcounter{revcomment}
    
    \edef\currentrevno{\the\value{reviewer}}
    \edef\currentcomno{\the\value{revcomment}}
    
    \edef\@currentlabel{\currentrevno.\currentcomno}
    \def\currlabel{rev:r\currentrevno c\currentcomno}
    \label{\currlabel}
    
	\reviewresponse@setprioritystyle{#1}
	\reviewresponse@setcompletedstyle{#2}

	\ifreviewresponselevel
		\ifreviewresponsecompleted
			\addcontentsline{toc}{subsection}{Comment \currentcomno{} - Done!}
		\else
			\addcontentsline{toc}{subsection}{Comment \currentcomno}
		\fi
	\else
		\addcontentsline{toc}{subsection}{Comment \currentcomno}
	\fi

	\ifreviewresponselevel
		\ifreviewresponsecompleted
			\def\reviewresponse@commenttitle{Reviewer \currentrevno{} Comment \currentcomno{} \hfill [\reviewresponse@prioritytitle{} - Done!]}
		\else
			\def\reviewresponse@commenttitle{Reviewer \currentrevno{} Comment \currentcomno{} \hfill [\reviewresponse@prioritytitle]}
		\fi
	\else
		\def\reviewresponse@commenttitle{Reviewer \currentrevno{} Comment \currentcomno{}}
	\fi

    \begin{tcolorbox}[
		adjusted title={\reviewresponse@commenttitle}, 
		fonttitle={\rmfamily\bfseries}, 
		colback={\reviewresponse@prioritybg}, 
		colframe={\reviewresponse@priorityframe},
        coltitle={white}, 
        bookmark={Reviewer \currentrevno{} Comment \currentcomno},
		#3
    ]
}{\end{tcolorbox}}

\newenvironment{howtoanswer}[1][{}]{%
	\textbf{How to answer:} #1\par
}{\vspace{1em plus 0.2em minus 1.5em}}

\newenvironment{revresponse}[1][]{%
	\textbf{Response:} #1\par
}{\vspace{4em plus 0.2em minus 1.5em}}

\newenvironment{concludingresponse}[1][]{%
	\begin{tcolorbox}[attach title to upper, breakable, title={Concluding Response #1}, after title={.\enskip}, fonttitle={\rmfamily\bfseries}, coltitle={colorconcludingresponsefg}, colback={colorconcludingresponsebg}, colframe={colorconcludingresponseframe},]
	}{\end{tcolorbox}}

\newenvironment{changes}{\begin{tcolorbox}[breakable,colback={colorchangebg}, colframe={colorchangeframe},enhanced jigsaw,]
	}{\end{tcolorbox}}

\newcommand{\printpartbibliography}[1]{\begin{refsegment}\nocite{#1}\printbibliography[heading=none,segment=\therefsegment]\end{refsegment}}


\providecommand*{\revcommentautorefname}{Comment}
\AddToHook{package/hyperref/after}{
	\hypersetup{
		bookmarksopen,
		pdfcreator={LaTeX with the reviewresponse package},
	}
}
\AddToHook{begindocument/before}{
	\IfPackageLoadedTF{hyperref}{
		\hypersetup{
			pdfauthor={\@author},
			pdftitle={Review Response -- \@title{}},
			pdfkeywords={\@journal{}, \@manuscript{}, review response}
		}
	}{}
}
\AddToHook{begindocument/end}{\def\sectionautorefname{Reviewer}}