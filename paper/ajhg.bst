% version = 1.03 of ajhg.bst 2008 Feb 20
% American Journal of Human Genetics
% BibTeX format
%
% Send bug reports to Amy L. Williams <amy@csail.mit.edu>.


%%%% BUGS
%%%% InCollection and InBook print pages before address/publisher!!!
%%%% If this bug bites you, email me -- it'll motivate me to fix it.
%%%% Better yet, see if you can fix it first.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change history
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2008 Feb 20 (1.03);
%     * Don't output comma in output.check for mid.sentence state
%	(Necessary to avoid printing comma between journal name and volume
%	 number for article entries.)
% 2007 Nov 12 (1.02):
%   Complete rewrite starting from the unsrt BibTeX style file.
%   Changes: (All marked with initials ALW)
%     * No more than 10 author names, followed by "et al."
%     * Author names listed as "Von Last, F., Jr."  (Von and Junior being
%	optional)
%     * Title as default in BibTeX style files where all bu the first letter of
%       a title (or the word following a colon) are capitalized.  To keep a word
%       in the case it is written in, surround it by braces.
%     * Only print the year for the date, no month; also print year in parens.
%     * Do not emphasize (italicize) book titles.
%     * Use pp. to designate page numbers (not pages and page) for all entry
%	types except articles.
%     * Don't include the "number" (as in volume(number), pages) for articles.
%     * Separate volume number from page numbers by a comma not a colon.
%     * Never includes chapter number.
%     * Book editors names are listed followed by ", eds." or ", ed." for
%	single.
%     * Year follows author list in all types.
%     * Do not emphasize (italicize) journal name.
%     * Do emphasize volume number.
%     * Address and Publisher should occur in that order, in parens, separated
%       by a colon -- this changed in multiple places.



% previous versions based on Chicago format: (version 1.02 was complete rewrite)
% 2000 July 5 (1.01) count of authors in bibliography (up to seven)
%                    count in label (up to two). style of authors.
%                    punctuation...
% 2000 May 18 (1.00) origin

% from
% BibTeX standard bibliography style `unsrt'
	% version 0.99a for BibTeX versions 0.99a or later, LaTeX version 2.09.
	% Copyright (C) 1985, all rights reserved.
	% Copying of this file is authorized only if either
	% (1) you make absolutely no changes to your copy, including name, or
	% (2) if you do make changes, you name it something other than
	% btxbst.doc, plain.bst, unsrt.bst, alpha.bst, and abbrv.bst.
	% This restriction helps ensure that all standard styles are identical.
	% The file btxbst.doc has the documentation for this style.

ENTRY
  { address
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
    month
    note
    number
    organization
    pages
    publisher
    school
    series
    title
    type
    volume
    year
  }
  {}
  { label }

INTEGERS { output.state before.all mid.sentence after.sentence after.block }

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}

STRINGS { s t }

FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    % Changed by ALW -- no comma here, extra punctuation is generally not
    % necessary.
    %{ ", " * write$ }
    { " " * write$ }
    { output.state after.block =
	{ add.period$ write$
	  newline$
	  "\newblock " write$
	}
	{ output.state before.all =
	    'write$
	    { add.period$ " " * write$ }
	  if$
	}
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

FUNCTION {output.bibitem}
{ newline$
  "\bibitem{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {fin.entry}
{ add.period$
  write$
  newline$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}

FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
	'skip$
	{ after.sentence 'output.state := }
      if$
    }
  if$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {new.block.checka}
{ empty$
    'skip$
    'new.block
  if$
}

FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}

FUNCTION {new.sentence.checka}
{ empty$
    'skip$
    'new.sentence
  if$
}

FUNCTION {new.sentence.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.sentence
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "{\em " swap$ * "}" * }
  if$
}

INTEGERS { nameptr namesleft numnames }

FUNCTION {format.names}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    % Default format changed by ALW:
%    { s nameptr "{ff~}{vv~}{ll}{, jj}" format.name$ 't :=
    %  Use Author format:  "Von {vv~} Last {ll}, F. {f.}, Junior"
    { s nameptr "{vv~}{ll}{, f.}{, jj}" format.name$ 't :=
      % This condition added by ALW: no more than 10 authors before et al
      nameptr #10 > % max number of authors in list
	{
	  ", et~al." *
	  % Note: we subtract 1 from this below, so this ends the loop
	  #1 'namesleft :=
	}
	{ % the following is from the original unsrt file:
	  nameptr #1 >
	    { namesleft #1 >
		{ ", " * t * }
		{ numnames #2 >
		    { "," * }
		    'skip$
		  if$
		  t "others" =
		    { " et~al." * }
		    { " and " * t * }
		  if$
		}
	      if$
	    }
	    't
	  if$
	}
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {format.authors}
{ author empty$
    { "" }
    { author format.names }
  if$
}

FUNCTION {format.editors}
{ editor empty$
    { "" }
    { editor format.names
      editor num.names$ #1 >
	% Changed by ALW -- list editors followed by ", eds." or ", ed.":
%	{ ", editors" * }
%	{ ", editor" * }
	{ ", eds." * }
	{ ", ed." * }
      if$
    }
  if$
}

FUNCTION {format.title}
{ title empty$
    { "" }
    { title "t" change.case$ }
  if$
}

FUNCTION {n.dashify}
{ 't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
	{ t #1 #2 substring$ "--" = not
	    { "--" *
	      t #2 global.max$ substring$ 't :=
	    }
	    {   { t #1 #1 substring$ "-" = }
		{ "-" *
		  t #2 global.max$ substring$ 't :=
		}
	      while$
	    }
	  if$
	}
	{ t #1 #1 substring$ *
	  t #2 global.max$ substring$ 't :=
	}
      if$
    }
  while$
}

FUNCTION {format.date}
{ year empty$
  { "" }
  % Changed by ALW -- we only print year and we print it in parens:
  { "(" year * ")" * }
%    { month empty$
%	{ "" }
%	{ "there's a month but no year in " cite$ * warning$
%	  month
%	}
%      if$
%    }
%    { month empty$
%	'year
%	{ month " " * year * }
%      if$
%    }
  if$
}

FUNCTION {format.btitle}
% Changed by ALW -- no emphasis
%{ title emphasize
{ title
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { "volume" volume tie.or.space.connect
      series empty$
	'skip$
	{ " of " * series emphasize * }
      if$
      "volume and number" number either.or.check
    }
  if$
}

FUNCTION {format.number.series}
{ volume empty$
    { number empty$
	{ series field.or.null }
	{ output.state mid.sentence =
	    { "number" }
	    { "Number" }
	  if$
	  number tie.or.space.connect
	  series empty$
	    { "there's a number but no series in " cite$ * warning$ }
	    { " in " * series * }
	  if$
	}
      if$
    }
    { "" }
  if$
}

FUNCTION {format.edition}
{ edition empty$
    { "" }
    { output.state mid.sentence =
	{ edition "l" change.case$ " edition" * }
	{ edition "t" change.case$ " edition" * }
      if$
    }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
	{ #1 'multiresult := }
	{ t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}

FUNCTION {format.pages}
{ pages empty$
    { "" }
    { pages multi.page.check
      % Changed by ALW -- pp. for single and multi-pages:
%	{ "pages" pages n.dashify tie.or.space.connect }
%	{ "page" pages tie.or.space.connect }
	{ "pp." pages n.dashify tie.or.space.connect }
	{ "pp." pages tie.or.space.connect }
      if$
    }
  if$
}

FUNCTION {format.vol.num.pages}
  % Changed by ALW -- don't print number (of volume); do emphasize volume:
{ volume field.or.null emphasize
%{ volume field.or.null
%  number empty$
%    'skip$
%    { "(" number * ")" * *
%      volume empty$
%	{ "there's a number but no volume in " cite$ * warning$ }
%	'skip$
%      if$
%    }
%  if$
  pages empty$
    'skip$
    { duplicate$ empty$
	{ pop$ format.pages }
	% Changed by ALW -- use comma followed by space instead of colon after
	% volume number to delineate from page numbers:
%	{ ":" * pages n.dashify * }
	{ ", " * pages n.dashify * }
      if$
    }
  if$
}

FUNCTION {format.chapter.pages}
% Changed by ALW -- don't include chapter number, just pages:
{ 'format.pages
%{ chapter empty$
%    'format.pages
%    { type empty$
%	{ "chapter" }
%	{ type "l" change.case$ }
%      if$
%      chapter tie.or.space.connect
%      pages empty$
%	'skip$
%	{ ", " * format.pages * }
%      if$
%    }
%  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { editor empty$
        % Changed by ALW -- do not emphasize book titles
%	{ "In " booktitle emphasize * }
%	{ "In " format.editors * ", " * booktitle emphasize * }
	{ "In " booktitle * }
	{ "In " format.editors * ", " * booktitle * }
      if$
    }
  if$
}

% Added by ALW -- print address and publisher in the following format
% (Address: Publisher)
FUNCTION {format.address.publisher}
{ address empty$
    { "(" publisher * ")" * }
    { "(" address * ": " * publisher * ")" * }
  if$
}

FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  month empty$ year empty$ note empty$
  and and and and and
    { "all relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type "t" change.case$
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { "Technical Report" }
    'type
  if$
  number empty$
    { "t" change.case$ }
    { number tie.or.space.connect }
  if$
}

FUNCTION {format.article.crossref}
{ key empty$
    { journal empty$
	{ "need key or journal for " cite$ * " to crossref " * crossref *
	  warning$
	  ""
	}
	{ "In {\em " journal * "\/}" * }
      if$
    }
    { "In " key * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.crossref.editor}
{ editor #1 "{vv~}{ll}" format.name$
  editor num.names$ duplicate$
  #2 >
    { pop$ " et~al." * }
    { #2 <
	'skip$
	{ editor #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
	    { " et~al." * }
	    { " and " * editor #2 "{vv~}{ll}" format.name$ * }
	  if$
	}
      if$
    }
  if$
}

FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      "In "
    }
    { "Volume" volume tie.or.space.connect
      " of " *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ series empty$
	    { "need editor, key, or series for " cite$ * " to crossref " *
	      crossref * warning$
	      "" *
	    }
	    { "{\em " * series * "\/}" * }
	  if$
	}
	{ key * }
      if$
    }
    { format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.incoll.inproc.crossref}
{ editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ booktitle empty$
	    { "need editor, key, or booktitle for " cite$ * " to crossref " *
	      crossref * warning$
	      ""
	    }
	    { "In {\em " booktitle * "\/}" * }
	  if$
	}
	{ "In " key * }
      if$
    }
    { "In " format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  new.block
  % Changed by ALW -- (next two lines added) put year after authors:
  format.date "year" output.check
  new.block
  format.title "title" output.check
  new.block
  crossref missing$
    % Changed by ALW -- do not emphasize journal name:
    %{ journal emphasize "journal" output.check
    { journal "journal" output.check
      format.vol.num.pages output
      % Changed by ALW -- year comes earlier:
      %format.date "year" output.check
    }
    { format.article.crossref output.nonnull
      format.pages output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {book}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check }
    { format.authors output.nonnull
      crossref missing$
	{ "author and editor" editor either.or.check }
	'skip$
      if$
    }
  if$
  new.block
  % Changed by ALW -- (next two lines added) put year after authors:
  format.date "year" output.check
  new.block
  format.btitle "title" output.check
  crossref missing$
    { format.bvolume output
      new.block
      format.number.series output
      new.sentence
      % Changed by ALW -- new format function for address and publisher:
      format.address.publisher "publisher" output.check
%      publisher "publisher" output.check
%      address output
    }
    { new.block
      format.book.crossref output.nonnull
    }
  if$
  format.edition output
  % Changed by ALW -- year comes earlier:
%  format.date "year" output.check
%  new.block
  note output
  fin.entry
}

FUNCTION {booklet}
{ output.bibitem
  format.authors output
  new.block
  % Changed by ALW -- (next two lines added) put year after authors:
  format.date output
  new.block
  format.title "title" output.check
  howpublished address new.block.checkb
  howpublished output
  address output
  % Changed by ALW -- year comes earlier:
%  format.date output
%  new.block
  note output
  fin.entry
}

FUNCTION {inbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check }
    { format.authors output.nonnull
      crossref missing$
	{ "author and editor" editor either.or.check }
	'skip$
      if$
    }
  if$
  new.block
  % Changed by ALW -- (next two lines added) put year after authors:
  format.date "year" output.check
  new.block
  format.btitle "title" output.check
  crossref missing$
    { format.bvolume output
      format.chapter.pages "chapter and pages" output.check
      new.block
      format.number.series output
      new.sentence
      % Changed by ALW -- new format function for address and publisher:
      format.address.publisher "publisher" output.check
%      publisher "publisher" output.check
%      address output
    }
    { format.chapter.pages "chapter and pages" output.check
      new.block
      format.book.crossref output.nonnull
    }
  if$
  format.edition output
  % Changed by ALW -- year comes earlier:
%  format.date "year" output.check
%  new.block
  note output
  fin.entry
}

FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  new.block
  % Changed by ALW -- (next two lines added) put year after authors:
  format.date "year" output.check
  new.block
  format.title "title" output.check
  new.block
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.bvolume output
      format.number.series output
      format.chapter.pages output
      new.sentence
      % Changed by ALW -- new format function for address and publisher:
      format.address.publisher "publisher" output.check
%      publisher "publisher" output.check
%      address output
      format.edition output
      % Changed by ALW -- year comes earlier:
%      format.date "year" output.check
    }
    { format.incoll.inproc.crossref output.nonnull
      format.chapter.pages output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  new.block
  % Changed by ALW -- (next two lines added) put year after authors:
  format.date "year" output.check
  new.block
  format.title "title" output.check
  new.block
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.bvolume output
      format.number.series output
      % Changed by ALW: address and publisher (if present) should preceed
      % pages
      publisher empty$
	{ organization output } % I have no idea what this is
	{ new.sentence
	  format.address.publisher output
	  organization output  % I have no idea what this is
	}
      if$
      % Note: the year used to come here, but now comes earlier
      format.pages output
%      format.pages output
%      address empty$
%	{ organization publisher new.sentence.checkb
%	  organization output
%	  publisher output
%	  format.date "year" output.check
%	}
%	{ address output.nonnull
%	  format.date "year" output.check
%	  new.sentence
%	  organization output
%	  publisher output
%	}
%      if$
    }
    { format.incoll.inproc.crossref output.nonnull
      format.pages output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  author empty$
    { organization empty$
	'skip$
	{ organization output.nonnull
	  address output
	}
      if$
    }
    { format.authors output.nonnull }
  if$
  new.block
  % Changed by ALW -- (next two lines added) put year after authors:
  format.date output
  new.block
  format.btitle "title" output.check
  author empty$
    { organization empty$
	{ address new.block.checka
	  address output
	}
	'skip$
      if$
    }
    { organization address new.block.checkb
      organization output
      address output
    }
  if$
  format.edition output
  % Changed by ALW -- year comes earlier:
%  format.date output
  new.block
  note output
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  new.block
  % Changed by ALW -- (next two lines added) put year after authors:
  format.date "year" output.check
  new.block
  format.title "title" output.check
  new.block
  "Master's thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  % Changed by ALW -- year comes earlier:
%  format.date "year" output.check
  new.block
  note output
  fin.entry
}

FUNCTION {misc}
{ output.bibitem
  format.authors output
  % Changed by ALW -- put year after authors:
  format.date output
  title howpublished new.block.checkb
  format.title output
  howpublished new.block.checka
  howpublished output
  % Changed by ALW -- year comes earlier:
%  format.date output
  new.block
  note output
  fin.entry
  empty.misc.check
}

FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  new.block
  % Changed by ALW -- (next two lines added) put year after authors:
  format.date "year" output.check
  new.block
  format.btitle "title" output.check
  new.block
  "PhD thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  % Changed by ALW -- year comes earlier:
%  format.date "year" output.check
  new.block
  note output
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  editor empty$
    { organization output }
    { format.editors output.nonnull }
  if$
  new.block
  % Changed by ALW -- (next two lines added) put year after authors:
  format.date "year" output.check
  new.block
  format.btitle "title" output.check
  format.bvolume output
  format.number.series output
  address empty$
    { editor empty$
	{ publisher new.sentence.checka }
	{ organization publisher new.sentence.checkb
	  organization output
	}
      if$
      publisher output
      % Changed by ALW -- year comes earlier:
      %format.date "year" output.check
    }
    { address output.nonnull
      % Changed by ALW -- year comes earlier:
      %format.date "year" output.check
      new.sentence
      editor empty$
	'skip$
	{ organization output }
      if$
      publisher output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  new.block
  % Changed by ALW -- (next two lines added) put year after authors:
  format.date "year" output.check
  new.block
  format.title "title" output.check
  new.block
  format.tr.number output.nonnull
  institution "institution" output.check
  address output
  % Changed by ALW -- year comes earlier:
%  format.date "year" output.check
  new.block
  note output
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.title "title" output.check
  new.block
  note "note" output.check
  format.date output
  fin.entry
}

FUNCTION {default.type} { misc }

MACRO {jan} {"January"}

MACRO {feb} {"February"}

MACRO {mar} {"March"}

MACRO {apr} {"April"}

MACRO {may} {"May"}

MACRO {jun} {"June"}

MACRO {jul} {"July"}

MACRO {aug} {"August"}

MACRO {sep} {"September"}

MACRO {oct} {"October"}

MACRO {nov} {"November"}

MACRO {dec} {"December"}

MACRO {acmcs} {"ACM Computing Surveys"}

MACRO {acta} {"Acta Informatica"}

MACRO {cacm} {"Communications of the ACM"}

MACRO {ibmjrd} {"IBM Journal of Research and Development"}

MACRO {ibmsj} {"IBM Systems Journal"}

MACRO {ieeese} {"IEEE Transactions on Software Engineering"}

MACRO {ieeetc} {"IEEE Transactions on Computers"}

MACRO {ieeetcad}
 {"IEEE Transactions on Computer-Aided Design of Integrated Circuits"}

MACRO {ipl} {"Information Processing Letters"}

MACRO {jacm} {"Journal of the ACM"}

MACRO {jcss} {"Journal of Computer and System Sciences"}

MACRO {scp} {"Science of Computer Programming"}

MACRO {sicomp} {"SIAM Journal on Computing"}

MACRO {tocs} {"ACM Transactions on Computer Systems"}

MACRO {tods} {"ACM Transactions on Database Systems"}

MACRO {tog} {"ACM Transactions on Graphics"}

MACRO {toms} {"ACM Transactions on Mathematical Software"}

MACRO {toois} {"ACM Transactions on Office Information Systems"}

MACRO {toplas} {"ACM Transactions on Programming Languages and Systems"}

MACRO {tcs} {"Theoretical Computer Science"}

READ

STRINGS { longest.label }

INTEGERS { number.label longest.label.width }

FUNCTION {initialize.longest.label}
{ "" 'longest.label :=
  #1 'number.label :=
  #0 'longest.label.width :=
}

FUNCTION {longest.label.pass}
{ number.label int.to.str$ 'label :=
  number.label #1 + 'number.label :=
  label width$ longest.label.width >
    { label 'longest.label :=
      label width$ 'longest.label.width :=
    }
    'skip$
  if$
}

EXECUTE {initialize.longest.label}

ITERATE {longest.label.pass}

FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{"  longest.label  * "}" * write$ newline$
}

EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}
