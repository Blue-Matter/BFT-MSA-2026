
function Figure(fig)
  local cap = fig.caption
  if not cap or not cap.long or #cap.long == 0 then
    return fig
  end

  -- cap.long is a list of blocks; stringify them
  local txt = pandoc.utils.stringify(cap.long)

  -- 1. Convert "Figure 1:" → "Figure 1."
  txt = txt:gsub("Figure%s+(%d+):", "Figure %1.")

  -- 2. Split into prefix "Figure 1. " and the rest
  local prefix, rest = txt:match("^(Figure%s+%d+%.%s*)(.*)$")
  if not prefix then
    return fig
  end

  -- 3. Build inlines: bold prefix + normal rest
  local inlines = {
    pandoc.Strong{ pandoc.Str(prefix) },
    pandoc.Str(rest)
  }

  -- 4. Wrap in a Plain block (because caption.long expects blocks)
  cap.long = { pandoc.Plain(inlines) }
  fig.caption = cap

  return fig
end

