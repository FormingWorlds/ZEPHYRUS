/* Render mermaid diagrams written as fenced code blocks.
 *
 * The site builder emits a ```mermaid fence as a highlighted code block, so
 * the diagram source reaches the page as text. This finds those blocks,
 * replaces each with the element mermaid renders into, and runs mermaid
 * once. Without the library (offline, or a blocked CDN) the page keeps the
 * code block and stays readable.
 */
(function () {
  var KEYWORDS = /^\s*(flowchart|graph|sequenceDiagram|classDiagram|stateDiagram(-v2)?|erDiagram|journey|gantt|pie|mindmap|timeline|quadrantChart|gitGraph)\b/;

  function collect() {
    // Both selectors can match inside one block, so key on the element
    // that gets replaced and keep the first source seen for it.
    var seen = new Map();
    document.querySelectorAll('pre code, div.highlight pre').forEach(function (node) {
      var text = node.textContent || '';
      if (!KEYWORDS.test(text)) {
        return;
      }
      var host = node.closest('div.highlight') || node.closest('pre');
      if (host && !seen.has(host)) {
        seen.set(host, text);
      }
    });
    return Array.from(seen, function (entry) {
      return { host: entry[0], source: entry[1] };
    });
  }

  function darkMode() {
    var scheme = document.body.getAttribute('data-md-color-scheme');
    if (scheme) {
      return scheme !== 'default';
    }
    return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
  }

  function render() {
    if (typeof window.mermaid === 'undefined') {
      return;
    }
    var blocks = collect();
    if (!blocks.length) {
      return;
    }
    blocks.forEach(function (block) {
      var target = document.createElement('pre');
      target.className = 'mermaid';
      target.textContent = block.source;
      block.host.replaceWith(target);
    });
    window.mermaid.initialize({
      startOnLoad: false,
      securityLevel: 'strict',
      theme: darkMode() ? 'dark' : 'default',
      flowchart: { htmlLabels: true, useMaxWidth: true },
    });
    window.mermaid.run({ querySelector: 'pre.mermaid' });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', render);
  } else {
    render();
  }
})();
