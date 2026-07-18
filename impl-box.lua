-- Render "### ... {.impl}" labels as a filled blue box in LaTeX (PDF).
-- HTML and EPUB keep the heading and style it via styles.css (h3.impl).
function Header(el)
  if FORMAT:match("latex") and el.classes:includes("impl") then
    local text = pandoc.utils.stringify(el)
    return pandoc.RawBlock("latex", "\\impllabel{" .. text .. "}")
  end
end
