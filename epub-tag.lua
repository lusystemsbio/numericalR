-- EPUB-only: make numbered equations render (with their number) in EPUB.
--
-- EPUB math is converted by Pandoc's texmath to MathML, whose parser chokes on
-- `\tag{N}` when a plain letter precedes it (e.g. `g - kX \tag{1}`): it fails the
-- whole equation and dumps the raw `$$...$$` / `\begin{equation}...` LaTeX into the
-- page. In the cases where it does parse, it silently drops the number. Both are
-- wrong for a book whose prose says "Equation (N)".
--
-- Fix: for EPUB output only, rewrite `\tag{N}` to `\qquad (N)`, which texmath
-- parses cleanly and renders as a visible "(N)" trailing the equation. This runs on
-- every carrier of equation TeX: `$$...$$` blocks arrive as Math elements, while
-- `\begin{equation}...\end{equation}` blocks arrive as raw tex (RawInline/RawBlock).
-- HTML (MathJax) and PDF (LaTeX amsmath) keep the real `\tag` and its flush-right
-- number, so this filter is a no-op for them.
if not FORMAT:match('epub') then
  return {}
end

local function retag(s)
  return (s:gsub('\\tag%*?%s*{([^}]*)}', '\\qquad (%1)'))
end

function Math(el)
  el.text = retag(el.text)
  return el
end

function RawInline(el)
  if el.format == 'tex' or el.format == 'latex' then
    el.text = retag(el.text)
    return el
  end
end

function RawBlock(el)
  if el.format == 'tex' or el.format == 'latex' then
    el.text = retag(el.text)
    return el
  end
end
