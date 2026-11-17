//! SVG file I/O and structured XML parsing.
//!
//! Provides [`write_svg_file`], [`read_svg_file`], and [`parse_svg`] backed by
//! `roxmltree`, returning a structured [`SvgDocument`] tree of [`SvgElement`]
//! nodes.

use std::path::Path;

/// A parsed SVG element node with its tag name, attributes, text content,
/// and child elements.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SvgElement {
    /// The XML tag name (for example `"rect"`, `"polyline"`, `"line"`, `"text"`).
    pub tag: String,
    /// Element attributes as `(name, value)` pairs in document order.
    pub attributes: Vec<(String, String)>,
    /// Direct non-whitespace text content of the element, if present.
    pub text: Option<String>,
    /// Direct child elements in document order.
    pub children: Vec<Self>,
}

impl SvgElement {
    /// Returns the value of the attribute named `name`, if present.
    #[must_use]
    pub fn attribute(
        &self,
        name: &str,
    ) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(k, _)| k == name)
            .map(|(_, v)| v.as_str())
    }
}

/// A parsed SVG document root with optional dimensions, `viewBox`, and child
/// elements.
#[derive(Debug, Clone, PartialEq)]
pub struct SvgDocument {
    /// Parsed numeric `width` attribute on the root `<svg>` element, if present.
    pub width: Option<f64>,
    /// Parsed numeric `height` attribute on the root `<svg>` element, if present.
    pub height: Option<f64>,
    /// Parsed `viewBox` (`[min_x, min_y, width, height]`), if present.
    pub view_box: Option<[f64; 4]>,
    /// Top-level child elements inside `<svg>`.
    pub elements: Vec<SvgElement>,
}

/// Parses a numeric SVG length (with optional `"px"` suffix) into a finite `f64`.
fn parse_svg_length(raw: &str) -> Result<f64, String> {
    let trimmed = raw.trim();
    let numeric = trimmed.strip_suffix("px").unwrap_or(trimmed).trim();
    let value: f64 = numeric
        .parse()
        .map_err(|e| format!("invalid SVG length {raw:?}: {e}"))?;
    if !value.is_finite() {
        return Err(format!("non-finite SVG length {raw:?}"));
    }
    Ok(value)
}

/// Converts a `roxmltree::Node` element into an [`SvgElement`] recursively.
#[allow(clippy::single_call_fn)]
fn element_from_node(node: roxmltree::Node<'_, '_>) -> SvgElement {
    let tag = node.tag_name().name().to_owned();
    let attributes = node
        .attributes()
        .map(|attr| (attr.name().to_owned(), attr.value().to_owned()))
        .collect();
    let text = node
        .text()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_owned);
    let children = node
        .children()
        .filter(roxmltree::Node::is_element)
        .map(element_from_node)
        .collect();
    SvgElement {
        tag,
        attributes,
        text,
        children,
    }
}

/// Parses an SVG string into a structured [`SvgDocument`].
///
/// # Errors
/// Fails when `svg` is not valid XML, the root element is not `<svg>`, or
/// `width`, `height`, or `viewBox` attributes contain invalid numbers.
pub fn parse_svg(svg: &str) -> Result<SvgDocument, String> {
    let doc = roxmltree::Document::parse(svg).map_err(|e| format!("invalid SVG XML: {e}"))?;
    let root = doc.root_element();
    if root.tag_name().name() != "svg" {
        return Err(format!(
            "expected <svg> root element, found <{}>",
            root.tag_name().name()
        ));
    }

    let width = root.attribute("width").map(parse_svg_length).transpose()?;
    let height = root.attribute("height").map(parse_svg_length).transpose()?;

    let view_box = match root.attribute("viewBox") {
        | Some(raw) => {
            let mut values = Vec::with_capacity(4);
            for part in raw
                .split(|c: char| c.is_ascii_whitespace() || c == ',')
                .filter(|s| !s.is_empty())
            {
                let v: f64 = part
                    .parse()
                    .map_err(|e| format!("invalid SVG viewBox number {part:?}: {e}"))?;
                if !v.is_finite() {
                    return Err(format!("non-finite SVG viewBox number {part:?}"));
                }
                values.push(v);
            }
            match values.as_slice() {
                | [x, y, w, h] => Some([*x, *y, *w, *h]),
                | _ => {
                    return Err(format!(
                        "expected 4 numbers in SVG viewBox, got {}",
                        values.len()
                    ));
                },
            }
        },
        | None => None,
    };

    let elements = root
        .children()
        .filter(roxmltree::Node::is_element)
        .map(element_from_node)
        .collect();

    Ok(SvgDocument {
        width,
        height,
        view_box,
        elements,
    })
}

/// Writes an SVG string to `path`.
///
/// # Errors
/// Fails when the file cannot be written.
pub fn write_svg_file(
    path: &Path,
    svg: &str,
) -> Result<(), String> {
    std::fs::write(path, svg).map_err(|e| format!("cannot write {}: {e}", path.display()))
}

/// Reads an SVG file at `path` and parses it into an [`SvgDocument`].
///
/// # Errors
/// Fails when the file cannot be read or is not valid SVG.
pub fn read_svg_file(path: &Path) -> Result<SvgDocument, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    parse_svg(&content)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_svg_extracts_dimensions_and_elements() -> Result<(), String> {
        let svg = r##"<svg xmlns="http://www.w3.org/2000/svg" width="800" height="600" viewBox="0 0 800 600">
<rect width="100%" height="100%" fill="#ffffff"/>
<text x="10" y="20">hello</text>
</svg>"##;
        let doc = parse_svg(svg)?;
        assert_eq!(doc.width, Some(800.0));
        assert_eq!(doc.height, Some(600.0));
        assert_eq!(doc.view_box, Some([0.0, 0.0, 800.0, 600.0]));
        assert_eq!(doc.elements.len(), 2);
        assert_eq!(doc.elements.first().map(|e| e.tag.as_str()), Some("rect"));
        assert_eq!(
            doc.elements.get(1).and_then(|e| e.text.as_deref()),
            Some("hello")
        );
        Ok(())
    }

    #[test]
    fn write_and_read_svg_file_round_trip() -> Result<(), String> {
        let path = std::env::temp_dir().join(format!("rssn_svg_test_{}.svg", std::process::id()));
        let svg = r#"<svg xmlns="http://www.w3.org/2000/svg" width="100" height="50"><line x1="0" y1="0" x2="100" y2="50"/></svg>"#;
        write_svg_file(&path, svg)?;
        let doc = read_svg_file(&path)?;
        let _r = std::fs::remove_file(&path);
        assert_eq!(doc.width, Some(100.0));
        assert_eq!(doc.height, Some(50.0));
        assert_eq!(doc.elements.len(), 1);
        Ok(())
    }

    #[test]
    fn parse_svg_rejects_non_svg_root() {
        assert!(parse_svg("<html/>").is_err());
    }
}
